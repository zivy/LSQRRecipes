#include <fstream>
#include <stdlib.h>
#include <vnl/vnl_cross.h>

#include "Frame.h"
#include "Point2D.h"
#include "Point3D.h"
#include "RandomNumberGenerator.h"
#include "PlanePhantomUSCalibrationParametersEstimator.h"

using namespace lsqrRecipes;

/**
 * Test using experimental data.
 * @param transformationsFileName File containing the rigid transformation data 
 *                                (T2i). Transformations are represented as a 
 *                                rotation matrix and translation. The expected
 *                                data representation is rotation matrix and 
 *                                translation vector. The order of the entries
 *                                in the file is (possibly on consecutive lines):
 *                                R00 R01 R02 t0 R10 R11 R12 t1 R20 R21 R22 t2
 *
 * @param imagePointsFileName File containing the 2D image point coordinates.
 *                            Each line in the file is:
 *                            px py
 */
int phantomTest(const std::string &transformationsFileName,
                const std::string &imagePointsFileName);

/**
 * Run the test for a planar phantom using simulated data. This is 
 * intended for regression testing.
 */
bool phantomTest();


bool loadTransformations(const std::string & transformationsFileName, 
                         std::vector<lsqrRecipes::Frame> &transformations);

bool loadImagePoints(const std::string &imagePointsFileName, 
                     size_t numberOfPoints, 
                     std::vector<lsqrRecipes::Point2D> &imagePoints);

void generatePhantomData(unsigned int numElements, double noiseStandardDeviation,
                         std::vector< lsqrRecipes::PlanePhantomUSCalibrationParametersEstimator::DataType > &data,
                         std::vector< lsqrRecipes::PlanePhantomUSCalibrationParametersEstimator::DataType > &cleanData,
                         std::vector<double> &usCalibrationParameters);

void reportResults(const std::string &title,                                     
                   const std::vector<double> &estimatedUSCalibrationParameters,
                   const std::vector< lsqrRecipes::PlanePhantomUSCalibrationParametersEstimator::DataType > &data);

bool reportAndCheckResults(const std::string &title,
                           const std::vector<double> &trueUSCalibrationParameters,
                           const std::vector<double> &estimatedUSCalibrationParameters,
                           const std::vector< lsqrRecipes::PlanePhantomUSCalibrationParametersEstimator::DataType > &data);


int main( int argc, char *argv[] )
{
  switch(argc) {
    case 1:  //test using simulated data      
      if(phantomTest())
        return EXIT_SUCCESS;
      else
        return EXIT_FAILURE;
    case 3:  //test with experimental data            
      return phantomTest(argv[1], argv[2]);
    default:
      std::cerr<<"Unexpected number of command line parameters.\n";
      std::cerr<<"Usage: \n";
      std::cerr<<"\t"<<argv[0]<<" transformationsFileName 2DPointsFileName\n";
  }
}
/*****************************************************************************/

int phantomTest(const std::string &transformationsFileName,
                const std::string &imagePointsFileName)


{
  size_t i, n;
  lsqrRecipes::PlanePhantomUSCalibrationParametersEstimator::DataType
    dataElement;
  std::vector< lsqrRecipes::PlanePhantomUSCalibrationParametersEstimator::DataType > data;  
 
  std::cout<<"Plane phantom, experimental data:\n";
  std::cout<<"---------------------------------\n";

    //load the calibration data
  std::vector<lsqrRecipes::Frame> transformations;
  std::vector<lsqrRecipes::Point2D> imagePoints;

  if(!loadTransformations(transformationsFileName, transformations)) {
    std::cerr<<"Failed to load transformations file.\n";
    return EXIT_FAILURE;
  }
  if(!loadImagePoints(imagePointsFileName, transformations.size(), imagePoints)) {
    std::cerr<<"Failed to load image points file.\n";
    return EXIT_FAILURE;
  }

  n = transformations.size();
  for(i=0; i<n; i++) {
    dataElement.T2 = transformations[i];
    dataElement.q = imagePoints[i];
    data.push_back(dataElement);
  }

  double maxDistanceBetweenPoints = 5.0;
  lsqrRecipes::PlanePhantomUSCalibrationParametersEstimator 
    usCalibration(maxDistanceBetweenPoints);

  std::vector<double> estimatedUSCalibrationParameters;
  
  usCalibration.setLeastSquaresType(lsqrRecipes::PlanePhantomUSCalibrationParametersEstimator::ANALYTIC);
  usCalibration.leastSquaresEstimate(data, estimatedUSCalibrationParameters);
  reportResults(std::string("Analytic least squares"), 
                estimatedUSCalibrationParameters,
                data);

  usCalibration.setLeastSquaresType(lsqrRecipes::PlanePhantomUSCalibrationParametersEstimator::ITERATIVE);
  usCalibration.leastSquaresEstimate(data, estimatedUSCalibrationParameters);
  reportResults(std::string("Iterative least squares"), 
                estimatedUSCalibrationParameters,
                data);
  return EXIT_SUCCESS;
}

/*****************************************************************************/
bool phantomTest()
{
  std::vector< lsqrRecipes::PlanePhantomUSCalibrationParametersEstimator::DataType > data, cleanData, minForEstimate;  
  std::vector<double> trueCalibrationParameters;
   bool res, ok = true;

  unsigned int NUM_ELEMENTS = 50;
  double STD_NOISE2D = 1.0;
  generatePhantomData(NUM_ELEMENTS, STD_NOISE2D, data, cleanData, trueCalibrationParameters);
 
  double maxDistanceBetweenPoints = 3.0;
  lsqrRecipes::PlanePhantomUSCalibrationParametersEstimator 
    usCalibration(maxDistanceBetweenPoints);

  std::vector<double> estimatedUSCalibrationParameters;

  std::cout<<"Wall phantom, simulated data:\n";
  std::cout<<"-----------------------------\n";
  reportResults(std::string("Known calibration parameters"), 
                trueCalibrationParameters,
                data);
          //estimate using minimal number of data elements, first 31 elements
          //of clean data
  std::vector< lsqrRecipes::PlanePhantomUSCalibrationParametersEstimator::DataType >::iterator it;
  it = cleanData.begin();
  for(size_t i=0; i<31; i++, it++);
  minForEstimate.insert(minForEstimate.begin(),cleanData.begin(), it);
  usCalibration.estimate(minForEstimate, estimatedUSCalibrationParameters);   
  ok = ok && reportAndCheckResults(std::string("Minimal number of elements, no noise"),
                                   trueCalibrationParameters,
                                   estimatedUSCalibrationParameters,
                                   data);  
      //check that the agree method works with clean data
  ok = ok && usCalibration.agree(estimatedUSCalibrationParameters, cleanData[0]);

         //analytic estimate, using noisy data
  usCalibration.setLeastSquaresType(lsqrRecipes::PlanePhantomUSCalibrationParametersEstimator::ANALYTIC);
  usCalibration.leastSquaresEstimate(data, estimatedUSCalibrationParameters);
            //don't use "ok = ok&&reportAndCheckResults()", if the compiler optimizes
            //for short circuit code then the function may not be called and the
            //information will not be printed
  res = reportAndCheckResults(std::string("Analytic least squares, noisy data"),
                              trueCalibrationParameters,
                              estimatedUSCalibrationParameters,
                              data);
  ok = ok && res;
          //iterative estimate, using noisy data
  usCalibration.setLeastSquaresType(lsqrRecipes::PlanePhantomUSCalibrationParametersEstimator::ITERATIVE);
  usCalibration.leastSquaresEstimate(data, estimatedUSCalibrationParameters);  
  res = reportAndCheckResults(std::string("Iterative least squares, noisy data"), 
                              trueCalibrationParameters,
                              estimatedUSCalibrationParameters,
                              data);
  ok = ok && res;
  return ok;
}

/*****************************************************************************/
bool loadTransformations(const std::string & transformationsFileName, 
                         std::vector<lsqrRecipes::Frame> &transformations)
{
  double R[3][3], t[3];
  lsqrRecipes::Frame f;

  transformations.clear();

  std::ifstream in;
  in.open(transformationsFileName.c_str());
  if(!in.is_open())
    return false;
  while(in>>R[0][0]>>R[0][1]>>R[0][2]>>t[0]>>R[1][0]>>R[1][1]>>R[1][2]>>t[1]>>R[2][0]>>R[2][1]>>R[2][2]>>t[2]) {
    f.setRotationMatrix(R);
    f.setTranslation(t);
    transformations.push_back(f);
  }
  in.close();

  if(!transformations.empty())
    return true;
  return false;
}

/*****************************************************************************/

bool loadImagePoints(const std::string &imagePointsFileName, 
                     size_t numberOfPoints, 
                     std::vector<lsqrRecipes::Point2D> &imagePoints)
{
  lsqrRecipes::Point2D p;

  imagePoints.clear();

  std::ifstream in;
  in.open(imagePointsFileName.c_str());
  if(!in.is_open())
    return false;

  while(in>>p[0]>>p[1] && imagePoints.size()<numberOfPoints) {
    imagePoints.push_back(p);
  }
  in.close();

  if(!imagePoints.empty() && imagePoints.size() == numberOfPoints)
    return true;
  imagePoints.clear();
  return false;
}

/*****************************************************************************/
void reportResults(const std::string &title,    
                   const std::vector<double> &estimatedUSCalibrationParameters,
                   const std::vector< lsqrRecipes::PlanePhantomUSCalibrationParametersEstimator::DataType > &data)
{
  std::vector<double> errors;
  double minErr, maxErr, meanErr;

  std::cout<<"**"<<title<<"**\n";  
  if(estimatedUSCalibrationParameters.size() == 0)
    std::cout<<"DEGENERATE CONFIGURATION\n\n\n";
  else {
    std::cout<<"omega1[y,x]:\n";
    std::cout<<"\t["<<estimatedUSCalibrationParameters[0]<<", "<<estimatedUSCalibrationParameters[1]<<"]\n";
    std::cout<<"t1[z]:\n";
    std::cout<<"\t["<<estimatedUSCalibrationParameters[2]<<"]\n";
    std::cout<<"t3[x,y,z]:\n";
    std::cout<<"\t["<<estimatedUSCalibrationParameters[3]<<", "<<estimatedUSCalibrationParameters[4];
    std::cout<<", "<<estimatedUSCalibrationParameters[5]<<"]\n";
    std::cout<<"omega3[z,y,x]:\n";
    std::cout<<"\t["<<estimatedUSCalibrationParameters[6]<<", "<<estimatedUSCalibrationParameters[7];
    std::cout<<", "<<estimatedUSCalibrationParameters[8]<<"]\n";
    std::cout<<"m[x,y]:\n";
    std::cout<<"\t["<<estimatedUSCalibrationParameters[9];
    std::cout<<", "<<estimatedUSCalibrationParameters[10]<<"]\n";

    lsqrRecipes::PlanePhantomUSCalibrationParametersEstimator::getDistanceStatistics(estimatedUSCalibrationParameters, data,
                                                                                     errors, minErr, maxErr, meanErr);
    double resNorm = 0.0;
    for(size_t z=0; z<errors.size(); z++)
    {
     resNorm+=(errors[z]*errors[z]);
    }

    std::cout<<"sum of squared errors: "<<resNorm<<"\n\n\n";
  }
}
/*****************************************************************************/
bool reportAndCheckResults(const std::string &title,
                           const std::vector<double> &trueUSCalibrationParameters,
                           const std::vector<double> &estimatedUSCalibrationParameters,
                           const std::vector< lsqrRecipes::PlanePhantomUSCalibrationParametersEstimator::DataType > &data)
{
  std::vector<double> errors;
  double minErr, maxErr, meanErr;
               //[mm]
  const double TRANSLATION_EPS = 3.0;
               //[degree]
  const double ANGULAR_EPS = 0.08726646259971647884618453842445;
             //[mm]/[pix]
  const double SCALE_EPS = 1.0; 

  std::cout<<"**"<<title<<"**\n";  
  if(estimatedUSCalibrationParameters.size() == 0) {
    std::cout<<"DEGENERATE CONFIGURATION\n\n\n";
    if(trueUSCalibrationParameters.size()!=0)
      return false;
    return true;
  }
  else {
    std::cout<<"omega1[y,x]:\n";
    std::cout<<"\t["<<estimatedUSCalibrationParameters[0];
    std::cout<<", "<<estimatedUSCalibrationParameters[1]<<"]\n";
    std::cout<<"t1[z]:\n";
    std::cout<<"\t["<<estimatedUSCalibrationParameters[2]<<"]\n";
    std::cout<<"t3[x,y,z]:\n";
    std::cout<<"\t["<<estimatedUSCalibrationParameters[3]<<", "<<estimatedUSCalibrationParameters[4];
    std::cout<<", "<<estimatedUSCalibrationParameters[5]<<"]\n";
    std::cout<<"omega3[z,y,x]:\n";
    std::cout<<"\t["<<estimatedUSCalibrationParameters[6]<<", "<<estimatedUSCalibrationParameters[7];
    std::cout<<", "<<estimatedUSCalibrationParameters[8]<<"]\n";
    std::cout<<"m[x,y]:\n";
    std::cout<<"\t["<<estimatedUSCalibrationParameters[9];
    std::cout<<", "<<estimatedUSCalibrationParameters[10]<<"]\n";

    lsqrRecipes::PlanePhantomUSCalibrationParametersEstimator::getDistanceStatistics(estimatedUSCalibrationParameters, data,
                                                                                     errors, minErr, maxErr, meanErr);
    double resNorm = 0.0;
    for(size_t z=0; z<errors.size(); z++)
    {
     resNorm+=(errors[z]*errors[z]);
    }
    std::cout<<"sum of squared errors: "<<resNorm<<"\n\n\n";

                //compare known and estimated values. checking the angular 
                //values requires extracting the two sets of Euler angles and 
                //then comparing to the known values
    vnl_matrix<double> R(3,3);
    vnl_vector<double> r1(3), r2(3), r3(3);
    double omega3_x1, omega3_x2, omega3_y1, omega3_y2, omega3_z1, omega3_z2, cy;
    double smallAngle = 0.008726535498373935;//0.5 degrees
    double halfPI = 1.5707963267948966192313216916398;

    r1(0) = estimatedUSCalibrationParameters[11]/(estimatedUSCalibrationParameters[9]*estimatedUSCalibrationParameters[38]);
    r1(1) = estimatedUSCalibrationParameters[12]/(estimatedUSCalibrationParameters[9]*estimatedUSCalibrationParameters[38]);
    r1(2) = estimatedUSCalibrationParameters[13]/(estimatedUSCalibrationParameters[9]*estimatedUSCalibrationParameters[38]);

    r2(0) = estimatedUSCalibrationParameters[20]/(estimatedUSCalibrationParameters[10]*estimatedUSCalibrationParameters[38]);
    r2(1) = estimatedUSCalibrationParameters[21]/(estimatedUSCalibrationParameters[10]*estimatedUSCalibrationParameters[38]);
    r2(2) = estimatedUSCalibrationParameters[22]/(estimatedUSCalibrationParameters[10]*estimatedUSCalibrationParameters[38]);
  
    r3 = vnl_cross_3d(r1,r2);
    R.set_column(0,r1);
    R.set_column(1,r2);    
    R.set_column(2,r3);
   
               //two options for omega_y depending on choice of +-1 for sqrt          
    omega3_y1 = atan2(-R(2,0), sqrt(R(0,0)*R(0,0) + R(1,0)*R(1,0)));
    omega3_y2 = atan2(-R(2,0), -sqrt(R(0,0)*R(0,0) + R(1,0)*R(1,0)));
             //omega_z and omega_x depend on the omega_y value
	  if(fabs(omega3_y1 - halfPI) > smallAngle && fabs(omega3_y1 + halfPI) > smallAngle) {
      cy = cos(omega3_y1);
      omega3_z1 = atan2(R(1,0)/cy, R(0,0)/cy);
      omega3_x1 = atan2(R(2,1)/cy, R(2,2)/cy);     
      cy = cos(omega3_y2);
      omega3_z2 = atan2(R(1,0)/cy, R(0,0)/cy);
      omega3_x2 = atan2(R(2,1)/cy, R(2,2)/cy);

    }     //gimbal lock, omega_y is approximatly plus/minus half PI, arbitrarily 
        //set omega_z to 0 and compute omega_x
    else {          
      omega3_z1 = omega3_z2 = 0;
      omega3_x1 = omega3_x2 = atan2(R(0,1), R(1,1));
    }

    bool angles3OK = (fabs(omega3_z1 - trueUSCalibrationParameters[6])<ANGULAR_EPS &&
                      fabs(omega3_y1 - trueUSCalibrationParameters[7])<ANGULAR_EPS &&
                      fabs(omega3_x1 - trueUSCalibrationParameters[8])<ANGULAR_EPS) ||
                     (fabs(omega3_z2 - trueUSCalibrationParameters[6])<ANGULAR_EPS &&
                      fabs(omega3_y2 - trueUSCalibrationParameters[7])<ANGULAR_EPS &&
                      fabs(omega3_x2 - trueUSCalibrationParameters[8])<ANGULAR_EPS);

                      //only check the parameters for T3, we don't care about T1
    return (fabs(estimatedUSCalibrationParameters[3] - trueUSCalibrationParameters[3])<TRANSLATION_EPS &&
            fabs(estimatedUSCalibrationParameters[4] - trueUSCalibrationParameters[4])<TRANSLATION_EPS &&
            fabs(estimatedUSCalibrationParameters[5] - trueUSCalibrationParameters[5])<TRANSLATION_EPS &&
            angles3OK &&
            fabs(estimatedUSCalibrationParameters[9] - trueUSCalibrationParameters[9])<SCALE_EPS &&
            fabs(estimatedUSCalibrationParameters[10] - trueUSCalibrationParameters[10])<SCALE_EPS);
  }
}
/*****************************************************************************/
void generatePhantomData(unsigned int numElements, double noiseStandardDeviation,
                         std::vector< lsqrRecipes::PlanePhantomUSCalibrationParametersEstimator::DataType > &data,
                         std::vector< lsqrRecipes::PlanePhantomUSCalibrationParametersEstimator::DataType > &cleanData,
                         std::vector<double> &usCalibrationParameters)
{
  vnl_matrix<double> T3(4,4), R2i(3,3), R1(3,3), T1inv(4,4);
  vnl_vector<double> t1(3), t1inv(3), t2i(3), piH(4), piTransformedH(4), piTransformed(3), qiH(4), qiTransformedH(4), qiTransformed(3);
  double omega3_x, omega3_y, omega3_z, omega2_x, omega2_y, omega2_z, omega1_x, omega1_y, omega1_z, m_x, m_y;
	double cx, cy, cz, sx, sy, sz;

  double translationBounds[] = {-100, 100};
  double rotationBounds[] = {0.0, 3.1415926535897932384626433832795};
  double imageSizes[] = {640, 480};

  RandomNumberGenerator random;


             //random T3 (US calibration transformation)
  m_x = 0.143;
  m_y = 0.139;

  omega3_x = random.uniform(rotationBounds[0], rotationBounds[1]);
  omega3_y = random.uniform(rotationBounds[0], rotationBounds[1]);
  omega3_z = random.uniform(rotationBounds[0], rotationBounds[1]);
	cx = cos(omega3_x);
	cy = cos(omega3_y);
	cz = cos(omega3_z);
	sx = sin(omega3_x);
	sy = sin(omega3_y);
	sz = sin(omega3_z);
  
	T3(0,0) = m_x*(cz*cy);   T3(0,1) = m_y*(cz*sy*sx - sz*cx);   T3(0,2) = cz*sy*cx+sz*sx;     T3(0,3) = random.uniform(translationBounds[0], translationBounds[1]);
  T3(1,0) = m_x*(sz*cy);   T3(1,1) = m_y*(sz*sy*sx + cz*cx);   T3(1,2) = sz*sy*cx - cz*sx;   T3(1,3) = random.uniform(translationBounds[0], translationBounds[1]);
  T3(2,0) = m_x*(-sy);     T3(2,1) = m_y*(cy*sx);              T3(2,2) = cy*cx;              T3(2,3) = random.uniform(translationBounds[0], translationBounds[1]);
  T3(3,0) = 0.0;           T3(3,1) = 0.0;                      T3(3,2) = 0.0;                T3(3,3) = 1.0;
               //random T1 (transformation from tracker/world coordinate system to phantom coordinate system)
  t1(0) = random.uniform(translationBounds[0], translationBounds[1]);
  t1(1) = random.uniform(translationBounds[0], translationBounds[1]);
  t1(2) = random.uniform(translationBounds[0], translationBounds[1]);
  omega1_x = random.uniform(rotationBounds[0], rotationBounds[1]);
  omega1_y = random.uniform(rotationBounds[0], rotationBounds[1]);
  omega1_z = random.uniform(rotationBounds[0], rotationBounds[1]);
	cx = cos(omega1_x);
	cy = cos(omega1_y);
	cz = cos(omega1_z);
	sx = sin(omega1_x);
	sy = sin(omega1_y);
	sz = sin(omega1_z);
	R1(0,0) = cz*cy;   R1(0,1) = cz*sy*sx - sz*cx;   R1(0,2) = cz*sy*cx+sz*sx;
  R1(1,0) = sz*cy;   R1(1,1) = sz*sy*sx + cz*cx;   R1(1,2) = sz*sy*cx - cz*sx;
  R1(2,0) = -sy;     R1(2,1) = cy*sx;              R1(2,2) = cy*cx;
  
           //first eleven parameters are the minimal parametrization, the rest are
           //dependent upon them, and are there to simplify other computations
  usCalibrationParameters.clear();
  usCalibrationParameters.push_back(omega1_y);
  usCalibrationParameters.push_back(omega1_x);
  usCalibrationParameters.push_back(t1(2));
  usCalibrationParameters.push_back(T3(0,3));
  usCalibrationParameters.push_back(T3(1,3));
  usCalibrationParameters.push_back(T3(2,3));
  usCalibrationParameters.push_back(omega3_z);
  usCalibrationParameters.push_back(omega3_y);
  usCalibrationParameters.push_back(omega3_x);
  usCalibrationParameters.push_back(m_x);
  usCalibrationParameters.push_back(m_y);
                       //carried along to simplify computations (don't want to 
                       //recompute these quantities multiple times)
  usCalibrationParameters.push_back(R1(2,0)*T3(0,0));
  usCalibrationParameters.push_back(R1(2,0)*T3(1,0));
  usCalibrationParameters.push_back(R1(2,0)*T3(2,0));
  usCalibrationParameters.push_back(R1(2,1)*T3(0,0));
  usCalibrationParameters.push_back(R1(2,1)*T3(1,0));
  usCalibrationParameters.push_back(R1(2,1)*T3(2,0));
  usCalibrationParameters.push_back(R1(2,2)*T3(0,0));
  usCalibrationParameters.push_back(R1(2,2)*T3(1,0));
  usCalibrationParameters.push_back(R1(2,2)*T3(2,0));
  usCalibrationParameters.push_back(R1(2,0)*T3(0,1));
  usCalibrationParameters.push_back(R1(2,0)*T3(1,1));
  usCalibrationParameters.push_back(R1(2,0)*T3(2,1));
  usCalibrationParameters.push_back(R1(2,1)*T3(0,1));
  usCalibrationParameters.push_back(R1(2,1)*T3(1,1));
  usCalibrationParameters.push_back(R1(2,1)*T3(2,1));
  usCalibrationParameters.push_back(R1(2,2)*T3(0,1));
  usCalibrationParameters.push_back(R1(2,2)*T3(1,1));
  usCalibrationParameters.push_back(R1(2,2)*T3(2,1));
  usCalibrationParameters.push_back(R1(2,0)*T3(0,3));
  usCalibrationParameters.push_back(R1(2,0)*T3(1,3));
  usCalibrationParameters.push_back(R1(2,0)*T3(2,3));
  usCalibrationParameters.push_back(R1(2,1)*T3(0,3));
  usCalibrationParameters.push_back(R1(2,1)*T3(1,3));
  usCalibrationParameters.push_back(R1(2,1)*T3(2,3));
  usCalibrationParameters.push_back(R1(2,2)*T3(0,3));
  usCalibrationParameters.push_back(R1(2,2)*T3(1,3));
  usCalibrationParameters.push_back(R1(2,2)*T3(2,3));
  usCalibrationParameters.push_back(R1(2,0));
  usCalibrationParameters.push_back(R1(2,1));
  usCalibrationParameters.push_back(R1(2,2));
  

  t1inv = - (R1.transpose() * t1);
	T1inv(0,0) = R1(0,0);   T1inv(0,1) = R1(1,0);   T1inv(0,2) = R1(2,0);   T1inv(0,3) = t1inv(0); 
  T1inv(1,0) = R1(0,1);   T1inv(1,1) = R1(1,1);   T1inv(1,2) = R1(2,1);   T1inv(1,3) = t1inv(1);
  T1inv(2,0) = R1(0,2);   T1inv(2,1) = R1(1,2);   T1inv(2,2) = R1(2,2);   T1inv(2,3) = t1inv(2);
  T1inv(3,0) = 0.0;       T1inv(3,1) = 0.0;       T1inv(3,2) = 0.0;       T1inv(3,3) = 1.0;


  lsqrRecipes::PlanePhantomUSCalibrationParametersEstimator::DataType d;
          //select a random point in the image and set the US probe's transformation
          //T2i accordingly, randomly generate a rotation and then the translation
          //is defined by this rotation+US calibration+target point in tracker
          //coordinate system, that is we have t2i(T3,t1,R2i).
  qiH[2] = 0.0;
  qiH[3] = 1.0;
  piH[2] = 0.0;
  piH[3] = 1.0;
  d.T2.setIdentity();
  for(unsigned int i=0; i<numElements; i++)
  {           //set the 2D image point
    qiH[0] = random.uniform(0,imageSizes[0]);
    qiH[1] = random.uniform(0,imageSizes[1]);
    d.q[0] = qiH[0];
    d.q[1] = qiH[1];
              //random point on the z=0 plane
    piH[0] = random.uniform(translationBounds[0], translationBounds[1]);
    piH[1] = random.uniform(translationBounds[0], translationBounds[1]);

           //get random rotation and dependent translation  
    omega2_x = random.uniform(rotationBounds[0], rotationBounds[1]);
    omega2_y = random.uniform(rotationBounds[0], rotationBounds[1]);
    omega2_z = random.uniform(rotationBounds[0], rotationBounds[1]);
	  cx = cos(omega2_x);
	  cy = cos(omega2_y);
	  cz = cos(omega2_z);
	  sx = sin(omega2_x);
	  sy = sin(omega2_y);
	  sz = sin(omega2_z);

    R2i(0,0) = cz*cy;   R2i(0,1) = cz*sy*sx - sz*cx;   R2i(0,2) = cz*sy*cx+sz*sx;
    R2i(1,0) = sz*cy;   R2i(1,1) = sz*sy*sx + cz*cx;   R2i(1,2) = sz*sy*cx - cz*sx;
    R2i(2,0) = -sy;     R2i(2,1) = cy*sx;              R2i(2,2) = cy*cx;
    
    qiTransformedH = T3*qiH;
    qiTransformed[0] = qiTransformedH[0];
    qiTransformed[1] = qiTransformedH[1];
    qiTransformed[2] = qiTransformedH[2];

    piTransformedH = T1inv*piH;
    piTransformed[0] = piTransformedH[0];
    piTransformed[1] = piTransformedH[1];
    piTransformed[2] = piTransformedH[2];
   
    t2i = piTransformed - R2i*qiTransformed;

    d.T2.setTranslation(t2i(0), t2i(1), t2i(2));
    d.T2.setRotationMatrix(R2i(0,0), R2i(0,1), R2i(0,2),
                           R2i(1,0), R2i(1,1), R2i(1,2),
                           R2i(2,0), R2i(2,1), R2i(2,2));

    cleanData.push_back(d);
           //add noise to the 2D points and save
    d.q[0] += random.normal(noiseStandardDeviation);
    d.q[1] += random.normal(noiseStandardDeviation);
    data.push_back(d);
  }
}
