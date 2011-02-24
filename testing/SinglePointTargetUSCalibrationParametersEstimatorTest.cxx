#include <fstream>
#include <stdlib.h>
#include "Frame.h"
#include "Point2D.h"
#include "Point3D.h"
#include "RandomNumberGenerator.h"
#include "SinglePointTargetUSCalibrationParametersEstimator.h"

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
 * Run the test for a single point phantom using simulated data. This is 
 * intended for regression testing.
 */
bool phantomTest();

/**
 * @param transformationsFileName File containing the rigid transformation data 
 *                                (T2i). Transformations are represented as a 
 *                                rotation matrix and translation. Each line in
 *                                the file is:
 *                                R00 R01 R02 t0 R10 R11 R12 t1 R20 R21 R22 t2
 *
 * @param imagePointsFileName File containing the 2D image point coordinates.
 *                            Each line in the file is:
 *                            px py
 */
int calibratedPointerTest(const std::string &transformationsFileName,
                          const std::string &imagePointsFileName,
                          const std::string &spatialPointsFileName);
bool calibratedPointerTest();

bool loadTransformations(const std::string & transformationsFileName, 
                         std::vector<lsqrRecipes::Frame> &transformations);

bool loadImagePoints(const std::string &imagePointsFileName, 
                     size_t numberOfPoints, 
                     std::vector<lsqrRecipes::Point2D> &imagePoints);

bool loadSpatialPoints(const std::string &spatialPointsFileName, 
                       size_t numberOfPoints, 
                       std::vector<lsqrRecipes::Point3D> &spatialPoints);

void generatePhantomData(unsigned int numElements, double noiseStandardDeviation,
                         std::vector< lsqrRecipes::SingleUnknownPointTargetUSCalibrationParametersEstimator::DataType > &data,
                         std::vector< lsqrRecipes::SingleUnknownPointTargetUSCalibrationParametersEstimator::DataType > &cleanData,
                         std::vector<double> &usCalibrationParameters);

void generateCalibratedPointerData(unsigned int numElements, double noiseStandardDeviation,
                                   std::vector< lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator::DataType > &data,
                                   std::vector< lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator::DataType > &cleanData,
                                   std::vector<double> &usCalibrationParameters);

void reportResultsSinglePointPhantom(const std::string &title,                                     
                                     const std::vector<double> &estimatedUSCalibrationParameters,
                                     const std::vector< lsqrRecipes::SingleUnknownPointTargetUSCalibrationParametersEstimator::DataType > &data);

bool reportAndCheckResultsSinglePointPhantom(const std::string &title,
                                             const std::vector<double> &trueUSCalibrationParameters,
                                             const std::vector<double> &estimatedUSCalibrationParameters,
                                             const std::vector< lsqrRecipes::SingleUnknownPointTargetUSCalibrationParametersEstimator::DataType > &data);

void reportResultsCalibratedPointer(const std::string &title,                                     
                                     const std::vector<double> &estimatedUSCalibrationParameters,
                                     const std::vector< lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator::DataType > &data);

bool reportAndCheckResultsCalibratedPointer(const std::string &title,
                                             const std::vector<double> &trueUSCalibrationParameters,
                                             const std::vector<double> &estimatedUSCalibrationParameters,
                                             const std::vector< lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator::DataType > &data);


int main( int argc, char *argv[] )
{
  switch(argc) {
    case 1:  //test using simulated data
      bool phantomTestResult, calibratedPointerTestResult; 
      phantomTestResult = phantomTest();     
      calibratedPointerTestResult = calibratedPointerTest();
      if(phantomTestResult && calibratedPointerTestResult)
        return EXIT_SUCCESS;
      else
        return EXIT_FAILURE;
    case 3:  //test with single target point phantom            
      return phantomTest(argv[1], argv[2]);
    case 4:  //test with calibrated pointer 
      return calibratedPointerTest(argv[1], argv[2], argv[3]);
    default:
      std::cerr<<"Unexpected number of command line parameters.\n";
      std::cerr<<"Usage: \n";
      std::cerr<<"\t"<<argv[0]<<" transformationsFileName 2DPointsFileName\n";
      std::cerr<<"\t\t\tor\n";
      std::cerr<<"\t"<<argv[0]<<" transformationsFileName 2DPointsFileName 3DPointsFileName\n";
  }
}
/*****************************************************************************/

int phantomTest(const std::string &transformationsFileName,
                const std::string &imagePointsFileName)


{
  size_t i, n;
  lsqrRecipes::SingleUnknownPointTargetUSCalibrationParametersEstimator::DataType
    dataElement;
  std::vector< lsqrRecipes::SingleUnknownPointTargetUSCalibrationParametersEstimator::DataType > data;  
 
  std::cout<<"Cross wire phantom, experimental data:\n";
  std::cout<<"--------------------------------------\n";

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
  lsqrRecipes::SingleUnknownPointTargetUSCalibrationParametersEstimator 
    usCalibration(maxDistanceBetweenPoints);

  std::vector<double> estimatedUSCalibrationParameters;
  
  usCalibration.setLeastSquaresType(lsqrRecipes::SingleUnknownPointTargetUSCalibrationParametersEstimator::ANALYTIC);
  usCalibration.leastSquaresEstimate(data, estimatedUSCalibrationParameters);
  reportResultsSinglePointPhantom(std::string("Analytic least squares"), 
                                  estimatedUSCalibrationParameters,
                                  data);

  usCalibration.setLeastSquaresType(lsqrRecipes::SingleUnknownPointTargetUSCalibrationParametersEstimator::ITERATIVE);
  usCalibration.leastSquaresEstimate(data, estimatedUSCalibrationParameters);
  reportResultsSinglePointPhantom(std::string("Iterative least squares"), 
                                  estimatedUSCalibrationParameters,
                                  data);
  return EXIT_SUCCESS;
}

/*****************************************************************************/
bool phantomTest()
{
  std::vector< lsqrRecipes::SingleUnknownPointTargetUSCalibrationParametersEstimator::DataType > data, cleanData, minForEstimate;  
  std::vector<double> trueCalibrationParameters;
  bool ok = true;

  unsigned int NUM_ELEMENTS = 50;
  double STD_NOISE2D = 1.0;
  generatePhantomData(NUM_ELEMENTS, STD_NOISE2D, data, cleanData, trueCalibrationParameters);
  
  double maxDistanceBetweenPoints = 3.0;
  lsqrRecipes::SingleUnknownPointTargetUSCalibrationParametersEstimator 
    usCalibration(maxDistanceBetweenPoints);

  std::vector<double> estimatedUSCalibrationParameters;

  std::cout<<"Cross wire phantom, simulated data:\n";
  std::cout<<"-----------------------------------\n";
  reportResultsSinglePointPhantom(std::string("Known calibration parameters"), 
                                  trueCalibrationParameters,
                                  data);
          //estimate using minimal number of data elements, first four elements
          //of clean data
  std::vector< lsqrRecipes::SingleUnknownPointTargetUSCalibrationParametersEstimator::DataType >::iterator it;
  it = cleanData.begin();
  for(size_t i=0; i<4; i++, it++);
  minForEstimate.insert(minForEstimate.begin(),cleanData.begin(), it);
  usCalibration.estimate(minForEstimate, estimatedUSCalibrationParameters);
  ok = ok && 
       reportAndCheckResultsSinglePointPhantom(std::string("Minimal number of elements, no noise"),
                                               trueCalibrationParameters,
                                               estimatedUSCalibrationParameters,
                                               data);
      //check that the agree method works with clean data
  ok = ok && usCalibration.agree(estimatedUSCalibrationParameters, cleanData[0]);

         //analytic estimate, using noisy data
  usCalibration.setLeastSquaresType(lsqrRecipes::SingleUnknownPointTargetUSCalibrationParametersEstimator::ANALYTIC);
  usCalibration.leastSquaresEstimate(data, estimatedUSCalibrationParameters);
  reportAndCheckResultsSinglePointPhantom(std::string("Analytic least squares, noisy data"),
                                  trueCalibrationParameters,
                                  estimatedUSCalibrationParameters,
                                  data);

          //iterative estimate, using noisy data
  usCalibration.setLeastSquaresType(lsqrRecipes::SingleUnknownPointTargetUSCalibrationParametersEstimator::ITERATIVE);
  usCalibration.leastSquaresEstimate(data, estimatedUSCalibrationParameters);
  ok = ok && 
       reportAndCheckResultsSinglePointPhantom(std::string("Iterative least squares, noisy data"), 
                                               trueCalibrationParameters,
                                               estimatedUSCalibrationParameters,
                                               data);
  return ok;
}

/*****************************************************************************/
int calibratedPointerTest(const std::string &transformationsFileName,
                          const std::string &imagePointsFileName,
                          const std::string &spatialPointsFileName)
{
  size_t i, n;
  std::vector<lsqrRecipes::Frame> transformations;
  std::vector<lsqrRecipes::Point2D> imagePoints;
  std::vector<lsqrRecipes::Point3D> spatialPoints;

  lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator::DataType
    dataElement;
  std::vector< lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator::DataType > data;  

  std::cout<<"Calibrated pointer, experimental data:\n";
  std::cout<<"--------------------------------------\n";

  if(!loadTransformations(transformationsFileName, transformations)) {
    std::cerr<<"Failed to load transformations file.\n";
    return EXIT_FAILURE;
  }
  if(!loadImagePoints(imagePointsFileName, transformations.size(), imagePoints)) {
    std::cerr<<"Failed to load image points file.\n";
    return EXIT_FAILURE;
  }
  if(!loadSpatialPoints(spatialPointsFileName, transformations.size(), spatialPoints)) {
    std::cerr<<"Failed to load image points file.\n";
    return EXIT_FAILURE;
  }

  n = transformations.size();
  for(i=0; i<n; i++) {
    dataElement.T2 = transformations[i];
    dataElement.q = imagePoints[i];
    dataElement.p = spatialPoints[i];
    data.push_back(dataElement);
  }

  double maxDistanceBetweenPoints = 5.0;
  lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator 
    usCalibration(maxDistanceBetweenPoints);

  std::vector<double> estimatedUSCalibrationParameters;
  
  usCalibration.setLeastSquaresType(lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator::ANALYTIC);
  usCalibration.leastSquaresEstimate(data, estimatedUSCalibrationParameters);
  reportResultsCalibratedPointer(std::string("Analytic least squares"), 
                                 estimatedUSCalibrationParameters,
                                 data);

  usCalibration.setLeastSquaresType(lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator::ITERATIVE);
  usCalibration.leastSquaresEstimate(data, estimatedUSCalibrationParameters);
  reportResultsCalibratedPointer(std::string("Iterative least squares"), 
                                  estimatedUSCalibrationParameters,
                                  data);
  return EXIT_SUCCESS;
}

/*****************************************************************************/
bool calibratedPointerTest()
{
  std::vector< lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator::DataType > data, cleanData, minForEstimate;  
  std::vector<double> trueCalibrationParameters;
  bool res, ok = true;

  unsigned int NUM_ELEMENTS = 50;
  double STD_NOISE2D = 1.0;
  generateCalibratedPointerData(NUM_ELEMENTS, STD_NOISE2D, data, cleanData, trueCalibrationParameters);
  
  std::cout<<"Calibrated pointer, simulated data:\n";
  std::cout<<"-----------------------------------\n";

  double maxDistanceBetweenPoints = 3.0;
  lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator 
    usCalibration(maxDistanceBetweenPoints);

  std::vector<double> estimatedUSCalibrationParameters;

  reportResultsCalibratedPointer(std::string("Known calibration parameters"), 
                                  trueCalibrationParameters,
                                  data);
          //estimate using minimal number of data elements, first four elements
          //of clean data
  std::vector< lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator::DataType >::iterator it;
  it = cleanData.begin();
  for(size_t i=0; i<3; i++, it++);
  minForEstimate.insert(minForEstimate.begin(),cleanData.begin(), it);
  usCalibration.estimate(minForEstimate, estimatedUSCalibrationParameters);
            //originally the line was ok = ok && reportAndCheckResultsCalibratedPointer(). As the
            //compiler performed short circuit evaluation when ok was already false the function
            //was not invoked and the information about failure was not reported
  res = reportAndCheckResultsCalibratedPointer(std::string("Minimal number of elements, no noise"),
                                               trueCalibrationParameters,
                                               estimatedUSCalibrationParameters,
                                               data);
  ok = ok && res; 

      //check that the agree method works with clean data
  res = usCalibration.agree(estimatedUSCalibrationParameters, cleanData[0]);
  ok = ok && res;

         //analytic estimate, using noisy data
  usCalibration.setLeastSquaresType(lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator::ANALYTIC);
  usCalibration.leastSquaresEstimate(data, estimatedUSCalibrationParameters);
  res = reportAndCheckResultsCalibratedPointer(std::string("Analytic least squares, noisy data"),
                                  trueCalibrationParameters,
                                  estimatedUSCalibrationParameters,
                                  data);
  ok = ok && res;
          //iterative estimate, using noisy data
  usCalibration.setLeastSquaresType(lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator::ITERATIVE);
  usCalibration.leastSquaresEstimate(data, estimatedUSCalibrationParameters);
  res = reportAndCheckResultsCalibratedPointer(std::string("Iterative least squares, noisy data"), 
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

bool loadSpatialPoints(const std::string &spatialPointsFileName, 
                       size_t numberOfPoints, 
                       std::vector<lsqrRecipes::Point3D> &spatialPoints)
{
  lsqrRecipes::Point3D p;

  spatialPoints.clear();

  std::ifstream in;
  in.open(spatialPointsFileName.c_str());
  if(!in.is_open())
    return false;

  while(in>>p[0]>>p[1]>>p[2] && spatialPoints.size()<numberOfPoints) {
    spatialPoints.push_back(p);
  }
  in.close();

  if(!spatialPoints.empty() && spatialPoints.size() == numberOfPoints)
    return true;
  spatialPoints.clear();
  return false;
}

/*****************************************************************************/
void reportResultsSinglePointPhantom(const std::string &title,    
                                     const std::vector<double> &estimatedUSCalibrationParameters,
                                     const std::vector< lsqrRecipes::SingleUnknownPointTargetUSCalibrationParametersEstimator::DataType > &data)
{
  std::vector<double> errors;
  double minErr, maxErr, meanErr;

  std::cout<<"**"<<title<<"**\n";  
  if(estimatedUSCalibrationParameters.size() == 0)
    std::cout<<"DEGENERATE CONFIGURATION\n\n\n";
  else {
    std::cout<<"t1[x,y,z]:\n";
    std::cout<<"\t["<<estimatedUSCalibrationParameters[0]<<", "<<estimatedUSCalibrationParameters[1];
    std::cout<<", "<<estimatedUSCalibrationParameters[2]<<"]\n";
    std::cout<<"t3[x,y,z]:\n";
    std::cout<<"\t["<<estimatedUSCalibrationParameters[3]<<", "<<estimatedUSCalibrationParameters[4];
    std::cout<<", "<<estimatedUSCalibrationParameters[5]<<"]\n";
    std::cout<<"omega[z,y,x]:\n";
    std::cout<<"\t["<<estimatedUSCalibrationParameters[6]<<", "<<estimatedUSCalibrationParameters[7];
    std::cout<<", "<<estimatedUSCalibrationParameters[8]<<"]\n";
    std::cout<<"m[x,y]:\n";
    std::cout<<"\t["<<estimatedUSCalibrationParameters[9];
    std::cout<<", "<<estimatedUSCalibrationParameters[10]<<"]\n";

    lsqrRecipes::SingleUnknownPointTargetUSCalibrationParametersEstimator::getDistanceStatistics(estimatedUSCalibrationParameters, data,
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
bool reportAndCheckResultsSinglePointPhantom(const std::string &title,
                                             const std::vector<double> &trueUSCalibrationParameters,
                                             const std::vector<double> &estimatedUSCalibrationParameters,
                                             const std::vector< lsqrRecipes::SingleUnknownPointTargetUSCalibrationParametersEstimator::DataType > &data)
{
  std::vector<double> errors;
  double minErr, maxErr, meanErr;
               //1[mm]
  const double TRANSLATION_EPS = 1.0;
               //1[degree]
  const double ANGULAR_EPS = 0.01745329251994329576923690768489;
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
    std::cout<<"t1[x,y,z]:\n";
    std::cout<<"\t["<<estimatedUSCalibrationParameters[0]<<", "<<estimatedUSCalibrationParameters[1];
    std::cout<<", "<<estimatedUSCalibrationParameters[2]<<"]\n";
    std::cout<<"t3[x,y,z]:\n";
    std::cout<<"\t["<<estimatedUSCalibrationParameters[3]<<", "<<estimatedUSCalibrationParameters[4];
    std::cout<<", "<<estimatedUSCalibrationParameters[5]<<"]\n";
    std::cout<<"omega[z,y,x]:\n";
    std::cout<<"\t["<<estimatedUSCalibrationParameters[6]<<", "<<estimatedUSCalibrationParameters[7];
    std::cout<<", "<<estimatedUSCalibrationParameters[8]<<"]\n";
    std::cout<<"m[x,y]:\n";
    std::cout<<"\t["<<estimatedUSCalibrationParameters[9];
    std::cout<<", "<<estimatedUSCalibrationParameters[10]<<"]\n";

    lsqrRecipes::SingleUnknownPointTargetUSCalibrationParametersEstimator::getDistanceStatistics(estimatedUSCalibrationParameters, data,
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
    double omega_x1, omega_x2, omega_y1, omega_y2, omega_z1, omega_z2, cy;
    double smallAngle = 0.008726535498373935;//0.5 degrees
    double halfPI = 1.5707963267948966192313216916398;

    R(0,0) = estimatedUSCalibrationParameters[11]/estimatedUSCalibrationParameters[9];
    R(1,0) = estimatedUSCalibrationParameters[12]/estimatedUSCalibrationParameters[9];
    R(2,0) = estimatedUSCalibrationParameters[13]/estimatedUSCalibrationParameters[9];
    R(0,1) = estimatedUSCalibrationParameters[14]/estimatedUSCalibrationParameters[10];
    R(1,1) = estimatedUSCalibrationParameters[15]/estimatedUSCalibrationParameters[10];
    R(2,1) = estimatedUSCalibrationParameters[16]/estimatedUSCalibrationParameters[10];
    R(0,2) = estimatedUSCalibrationParameters[17];
    R(1,2) = estimatedUSCalibrationParameters[18];
    R(2,2) = estimatedUSCalibrationParameters[19];
 
               //two options for omega_y depending on choice of +-1 for sqrt          
    omega_y1 = atan2(-R(2,0), sqrt(R(0,0)*R(0,0) + R(1,0)*R(1,0)));
    omega_y2 = atan2(-R(2,0), -sqrt(R(0,0)*R(0,0) + R(1,0)*R(1,0)));
             //omega_z and omega_x depend on the omega_y value
	  if(fabs(omega_y1 - halfPI) > smallAngle && fabs(omega_y1 + halfPI) > smallAngle) {
      cy = cos(omega_y1);
      omega_z1 = atan2(R(1,0)/cy, R(0,0)/cy);
      omega_x1 = atan2(R(2,1)/cy, R(2,2)/cy);
      cy = cos(omega_y2);
      omega_z2 = atan2(R(1,0)/cy, R(0,0)/cy);
      omega_x2 = atan2(R(2,1)/cy, R(2,2)/cy);

    }     //gimbal lock, omega_y is approximatly plus/minus half PI, arbitrarily 
        //set omega_z to 0 and compute omega_x
    else {          
      omega_z1 = omega_z2 = 0;
      omega_x1 = omega_x2 = atan2(R(0,1), R(1,1));
    }
             //check if one of the two sets of angles is close enough 
    bool anglesOK = (fabs(omega_z1 - trueUSCalibrationParameters[6])<ANGULAR_EPS &&
                   fabs(omega_y1 - trueUSCalibrationParameters[7])<ANGULAR_EPS &&
                   fabs(omega_x1 - trueUSCalibrationParameters[8])<ANGULAR_EPS) ||
                   (fabs(omega_z2 - trueUSCalibrationParameters[6])<ANGULAR_EPS &&
                   fabs(omega_y2 - trueUSCalibrationParameters[7])<ANGULAR_EPS &&
                   fabs(omega_x2 - trueUSCalibrationParameters[8])<ANGULAR_EPS);

                 //only check the parameters for T3, we don't care about T1
    return (fabs(estimatedUSCalibrationParameters[3] - trueUSCalibrationParameters[3])<TRANSLATION_EPS &&
            fabs(estimatedUSCalibrationParameters[4] - trueUSCalibrationParameters[4])<TRANSLATION_EPS &&
            fabs(estimatedUSCalibrationParameters[5] - trueUSCalibrationParameters[5])<TRANSLATION_EPS &&
            anglesOK &&
            fabs(estimatedUSCalibrationParameters[9] - trueUSCalibrationParameters[9])<SCALE_EPS &&
            fabs(estimatedUSCalibrationParameters[10] - trueUSCalibrationParameters[10])<SCALE_EPS);
  }
}
/*****************************************************************************/
void generatePhantomData(unsigned int numElements, double noiseStandardDeviation,
                         std::vector< lsqrRecipes::SingleUnknownPointTargetUSCalibrationParametersEstimator::DataType > &data,
                         std::vector< lsqrRecipes::SingleUnknownPointTargetUSCalibrationParametersEstimator::DataType > &cleanData,
                         std::vector<double> &usCalibrationParameters)
{
  vnl_matrix<double> T3(4,4), R2i(3,3);
  vnl_vector<double> t1(3), t2i(3), qiH(4), qiTransformedH(4), qiTransformed(3);
  double omega_x, omega_y, omega_z, m_x, m_y;
	double cx, cy, cz, sx, sy, sz;

  double translationBounds[] = {-100, 100};
  double rotationBounds[] = {0.0, 3.1415926535897932384626433832795};
  double imageSizes[] = {640, 480};

  RandomNumberGenerator random;


             //random T3 (US calibration transformation), and t1 (target point 
             //coordinates in tracker/world coordinate system)
  m_x = 0.143;
  m_y = 0.139;

  omega_x = random.uniform(rotationBounds[0], rotationBounds[1]);
  omega_y = random.uniform(rotationBounds[0], rotationBounds[1]);
  omega_z = random.uniform(rotationBounds[0], rotationBounds[1]);
	cx = cos(omega_x);
	cy = cos(omega_y);
	cz = cos(omega_z);
	sx = sin(omega_x);
	sy = sin(omega_y);
	sz = sin(omega_z);
  
	T3(0,0) = m_x*(cz*cy);   T3(0,1) = m_y*(cz*sy*sx - sz*cx);   T3(0,2) = cz*sy*cx+sz*sx;     T3(0,3) = random.uniform(translationBounds[0], translationBounds[1]);
  T3(1,0) = m_x*(sz*cy);   T3(1,1) = m_y*(sz*sy*sx + cz*cx);   T3(1,2) = sz*sy*cx - cz*sx;   T3(1,3) = random.uniform(translationBounds[0], translationBounds[1]);
  T3(2,0) = m_x*(-sy);     T3(2,1) = m_y*(cy*sx);              T3(2,2) = cy*cx;              T3(2,3) = random.uniform(translationBounds[0], translationBounds[1]);
  T3(3,0) = 0.0;           T3(3,1) = 0.0;                      T3(3,2) = 0.0;                T3(3,3) = 1.0;

  t1(0) = random.uniform(translationBounds[0], translationBounds[1]);
  t1(1) = random.uniform(translationBounds[0], translationBounds[1]);
  t1(2) = random.uniform(translationBounds[0], translationBounds[1]);

  usCalibrationParameters.clear();
  usCalibrationParameters.push_back(t1(0));
  usCalibrationParameters.push_back(t1(1));
  usCalibrationParameters.push_back(t1(2));
  usCalibrationParameters.push_back(T3(0,3));
  usCalibrationParameters.push_back(T3(1,3));
  usCalibrationParameters.push_back(T3(2,3));
  usCalibrationParameters.push_back(omega_z);
  usCalibrationParameters.push_back(omega_y);
  usCalibrationParameters.push_back(omega_x);
  usCalibrationParameters.push_back(m_x);
  usCalibrationParameters.push_back(m_y);
  usCalibrationParameters.push_back(T3(0,0));
  usCalibrationParameters.push_back(T3(1,0));
  usCalibrationParameters.push_back(T3(2,0));
  usCalibrationParameters.push_back(T3(0,1));
  usCalibrationParameters.push_back(T3(1,1));
  usCalibrationParameters.push_back(T3(2,1));
  usCalibrationParameters.push_back(T3(0,2));
  usCalibrationParameters.push_back(T3(1,2));
  usCalibrationParameters.push_back(T3(2,2));
  
  lsqrRecipes::SingleUnknownPointTargetUSCalibrationParametersEstimator::DataType d;
          //select a random point in the image and set the US probe's transformation
          //T2i accordingly, randomly generate a rotation and then the translation
          //is defined by this rotation+US calibration+target point in tracker
          //coordinate system, that is we have t2i(T3,t1,R2i).
  qiH[2] = 0.0;
  qiH[3] = 1.0;
  d.T2.setIdentity();
  for(unsigned int i=0; i<numElements; i++)
  {           //set the 2D image point
    qiH[0] = random.uniform(0,imageSizes[0]);
    qiH[1] = random.uniform(0,imageSizes[1]);
    d.q[0] = qiH[0];
    d.q[1] = qiH[1];
         
            //get random rotation and dependent translation  
    omega_x = random.uniform(rotationBounds[0], rotationBounds[1]);
    omega_y = random.uniform(rotationBounds[0], rotationBounds[1]);
    omega_z = random.uniform(rotationBounds[0], rotationBounds[1]);
	  cx = cos(omega_x);
	  cy = cos(omega_y);
	  cz = cos(omega_z);
	  sx = sin(omega_x);
	  sy = sin(omega_y);
	  sz = sin(omega_z);

    R2i(0,0) = cz*cy;   R2i(0,1) = cz*sy*sx - sz*cx;   R2i(0,2) = cz*sy*cx+sz*sx;
    R2i(1,0) = sz*cy;   R2i(1,1) = sz*sy*sx + cz*cx;   R2i(1,2) = sz*sy*cx - cz*sx;
    R2i(2,0) = -sy;     R2i(2,1) = cy*sx;              R2i(2,2) = cy*cx;
    
    qiTransformedH = T3*qiH;
    qiTransformed[0] = qiTransformedH[0];
    qiTransformed[1] = qiTransformedH[1];
    qiTransformed[2] = qiTransformedH[2];
    t2i = t1 - R2i*qiTransformed;

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
/*****************************************************************************/
void generateCalibratedPointerData(unsigned int numElements, double noiseStandardDeviation,
                                   std::vector< lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator::DataType > &data,
                                   std::vector< lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator::DataType > &cleanData,
                                   std::vector<double> &usCalibrationParameters)
{
  vnl_matrix<double> T3(4,4), T2i(4,4);
  vnl_vector<double> qiH(4), pi(4);
  double omega_x, omega_y, omega_z, m_x, m_y;
	double cx, cy, cz, sx, sy, sz;

  double translationBounds[] = {-100, 100};
  double rotationBounds[] = {0.0, 3.1415926535897932384626433832795};
  double imageSizes[] = {640, 480};

  RandomNumberGenerator random;


             //random T3 (US calibration transformation)
  m_x = 0.143;
  m_y = 0.139;

  omega_x = random.uniform(rotationBounds[0], rotationBounds[1]);
  omega_y = random.uniform(rotationBounds[0], rotationBounds[1]);
  omega_z = random.uniform(rotationBounds[0], rotationBounds[1]);
	cx = cos(omega_x);
	cy = cos(omega_y);
	cz = cos(omega_z);
	sx = sin(omega_x);
	sy = sin(omega_y);
	sz = sin(omega_z);
  
	T3(0,0) = m_x*(cz*cy);   T3(0,1) = m_y*(cz*sy*sx - sz*cx);   T3(0,2) = cz*sy*cx+sz*sx;     T3(0,3) = random.uniform(translationBounds[0], translationBounds[1]);
  T3(1,0) = m_x*(sz*cy);   T3(1,1) = m_y*(sz*sy*sx + cz*cx);   T3(1,2) = sz*sy*cx - cz*sx;   T3(1,3) = random.uniform(translationBounds[0], translationBounds[1]);
  T3(2,0) = m_x*(-sy);     T3(2,1) = m_y*(cy*sx);              T3(2,2) = cy*cx;              T3(2,3) = random.uniform(translationBounds[0], translationBounds[1]);
  T3(3,0) = 0.0;           T3(3,1) = 0.0;                      T3(3,2) = 0.0;                T3(3,3) = 1.0;

  usCalibrationParameters.clear();
  usCalibrationParameters.push_back(T3(0,3));
  usCalibrationParameters.push_back(T3(1,3));
  usCalibrationParameters.push_back(T3(2,3));
  usCalibrationParameters.push_back(omega_z);
  usCalibrationParameters.push_back(omega_y);
  usCalibrationParameters.push_back(omega_x);
  usCalibrationParameters.push_back(m_x);
  usCalibrationParameters.push_back(m_y);
  usCalibrationParameters.push_back(T3(0,0));
  usCalibrationParameters.push_back(T3(1,0));
  usCalibrationParameters.push_back(T3(2,0));
  usCalibrationParameters.push_back(T3(0,1));
  usCalibrationParameters.push_back(T3(1,1));
  usCalibrationParameters.push_back(T3(2,1));
  usCalibrationParameters.push_back(T3(0,2));
  usCalibrationParameters.push_back(T3(1,2));
  usCalibrationParameters.push_back(T3(2,2));
  
  lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator::DataType d;
          //select a random point in the image, a random US probe transformation
          //T2i, and the pointer tip location in the tracker coordinate system
          //that depends on these quantities pi(T3,T2i,qi).
  qiH[2] = 0.0;
  qiH[3] = 1.0;
  d.T2.setIdentity();
  for(unsigned int i=0; i<numElements; i++)
  {           //set the 2D image point
    qiH[0] = random.uniform(0,imageSizes[0]);
    qiH[1] = random.uniform(0,imageSizes[1]);
    d.q[0] = qiH[0];
    d.q[1] = qiH[1];
         
            //get random rotation and dependent 3D point  
    omega_x = random.uniform(rotationBounds[0], rotationBounds[1]);
    omega_y = random.uniform(rotationBounds[0], rotationBounds[1]);
    omega_z = random.uniform(rotationBounds[0], rotationBounds[1]);
	  cx = cos(omega_x);
	  cy = cos(omega_y);
	  cz = cos(omega_z);
	  sx = sin(omega_x);
	  sy = sin(omega_y);
	  sz = sin(omega_z);

	  T2i(0,0) = m_x*(cz*cy);   T2i(0,1) = m_y*(cz*sy*sx - sz*cx);   T2i(0,2) = cz*sy*cx+sz*sx;     T2i(0,3) = random.uniform(translationBounds[0], translationBounds[1]);
    T2i(1,0) = m_x*(sz*cy);   T2i(1,1) = m_y*(sz*sy*sx + cz*cx);   T2i(1,2) = sz*sy*cx - cz*sx;   T2i(1,3) = random.uniform(translationBounds[0], translationBounds[1]);
    T2i(2,0) = m_x*(-sy);     T2i(2,1) = m_y*(cy*sx);              T2i(2,2) = cy*cx;              T2i(2,3) = random.uniform(translationBounds[0], translationBounds[1]);
    T2i(3,0) = 0.0;           T2i(3,1) = 0.0;                      T2i(3,2) = 0.0;                T2i(3,3) = 1.0;

    d.T2.setTranslation(T2i(0,3), T2i(1,3), T2i(2,3));
    d.T2.setRotationMatrix(T2i(0,0), T2i(0,1), T2i(0,2),
                           T2i(1,0), T2i(1,1), T2i(1,2),
                           T2i(2,0), T2i(2,1), T2i(2,2));

    pi = T2i*T3*qiH;
    d.p[0] = pi[0];
    d.p[1] = pi[1];
    d.p[2] = pi[2];
    cleanData.push_back(d);
           //add noise to the 2D points and save
    d.q[0] += random.normal(noiseStandardDeviation);
    d.q[1] += random.normal(noiseStandardDeviation);
    data.push_back(d);
  }
}
/*****************************************************************************/
void reportResultsCalibratedPointer(const std::string &title,                                     
                                     const std::vector<double> &estimatedUSCalibrationParameters,
                                     const std::vector< lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator::DataType > &data)
{
  std::vector<double> errors;
  double minErr, maxErr, meanErr;

  std::cout<<"**"<<title<<"**\n";  
  if(estimatedUSCalibrationParameters.size() == 0)
    std::cout<<"DEGENERATE CONFIGURATION\n\n\n";
  else {
    std::cout<<"t3[x,y,z]:\n";
    std::cout<<"\t["<<estimatedUSCalibrationParameters[0]<<", "<<estimatedUSCalibrationParameters[1];
    std::cout<<", "<<estimatedUSCalibrationParameters[2]<<"]\n";
    std::cout<<"omega[z,y,x]:\n";
    std::cout<<"\t["<<estimatedUSCalibrationParameters[3]<<", "<<estimatedUSCalibrationParameters[4];
    std::cout<<", "<<estimatedUSCalibrationParameters[5]<<"]\n";
    std::cout<<"m[x,y]:\n";
    std::cout<<"\t["<<estimatedUSCalibrationParameters[6];
    std::cout<<", "<<estimatedUSCalibrationParameters[7]<<"]\n";

    lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator::getDistanceStatistics(estimatedUSCalibrationParameters, data,
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
bool reportAndCheckResultsCalibratedPointer(const std::string &title,
                                             const std::vector<double> &trueUSCalibrationParameters,
                                             const std::vector<double> &estimatedUSCalibrationParameters,
                                             const std::vector< lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator::DataType > &data)
{
  std::vector<double> errors;
  double minErr, maxErr, meanErr;
               //1[mm]
  const double TRANSLATION_EPS = 1.0;
               //1[degree]
  const double ANGULAR_EPS = 0.01745329251994329576923690768489;
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
    std::cout<<"t3[x,y,z]:\n";
    std::cout<<"\t["<<estimatedUSCalibrationParameters[0]<<", "<<estimatedUSCalibrationParameters[1];
    std::cout<<", "<<estimatedUSCalibrationParameters[2]<<"]\n";
    std::cout<<"omega[z,y,x]:\n";
    std::cout<<"\t["<<estimatedUSCalibrationParameters[3]<<", "<<estimatedUSCalibrationParameters[4];
    std::cout<<", "<<estimatedUSCalibrationParameters[5]<<"]\n";
    std::cout<<"m[x,y]:\n";
    std::cout<<"\t["<<estimatedUSCalibrationParameters[6];
    std::cout<<", "<<estimatedUSCalibrationParameters[7]<<"]\n";

    lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator::getDistanceStatistics(estimatedUSCalibrationParameters, data,
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
    double omega_x1, omega_x2, omega_y1, omega_y2, omega_z1, omega_z2, cy;
    double smallAngle = 0.008726535498373935;//0.5 degrees
    double halfPI = 1.5707963267948966192313216916398;

    R(0,0) = estimatedUSCalibrationParameters[8]/estimatedUSCalibrationParameters[6];
    R(1,0) = estimatedUSCalibrationParameters[9]/estimatedUSCalibrationParameters[6];
    R(2,0) = estimatedUSCalibrationParameters[10]/estimatedUSCalibrationParameters[6];
    R(0,1) = estimatedUSCalibrationParameters[11]/estimatedUSCalibrationParameters[7];
    R(1,1) = estimatedUSCalibrationParameters[12]/estimatedUSCalibrationParameters[7];
    R(2,1) = estimatedUSCalibrationParameters[13]/estimatedUSCalibrationParameters[7];
    R(0,2) = estimatedUSCalibrationParameters[14];
    R(1,2) = estimatedUSCalibrationParameters[15];
    R(2,2) = estimatedUSCalibrationParameters[16];
 
               //two options for omega_y depending on choice of +-1 for sqrt          
    omega_y1 = atan2(-R(2,0), sqrt(R(0,0)*R(0,0) + R(1,0)*R(1,0)));
    omega_y2 = atan2(-R(2,0), -sqrt(R(0,0)*R(0,0) + R(1,0)*R(1,0)));
             //omega_z and omega_x depend on the omega_y value
	  if(fabs(omega_y1 - halfPI) > smallAngle && fabs(omega_y1 + halfPI) > smallAngle) {
      cy = cos(omega_y1);
      omega_z1 = atan2(R(1,0)/cy, R(0,0)/cy);
      omega_x1 = atan2(R(2,1)/cy, R(2,2)/cy);
      cy = cos(omega_y2);
      omega_z2 = atan2(R(1,0)/cy, R(0,0)/cy);
      omega_x2 = atan2(R(2,1)/cy, R(2,2)/cy);

    }     //gimbal lock, omega_y is approximatly plus/minus half PI, arbitrarily 
        //set omega_z to 0 and compute omega_x
    else {          
      omega_z1 = omega_z2 = 0;
      omega_x1 = omega_x2 = atan2(R(0,1), R(1,1));
    }
             //check if one of the two sets of angles is close enough 
    bool anglesOK = (fabs(omega_z1 - trueUSCalibrationParameters[3])<ANGULAR_EPS &&
                   fabs(omega_y1 - trueUSCalibrationParameters[4])<ANGULAR_EPS &&
                   fabs(omega_x1 - trueUSCalibrationParameters[5])<ANGULAR_EPS) ||
                   (fabs(omega_z2 - trueUSCalibrationParameters[3])<ANGULAR_EPS &&
                   fabs(omega_y2 - trueUSCalibrationParameters[4])<ANGULAR_EPS &&
                   fabs(omega_x2 - trueUSCalibrationParameters[5])<ANGULAR_EPS);

    return (fabs(estimatedUSCalibrationParameters[0] - trueUSCalibrationParameters[0])<TRANSLATION_EPS &&
            fabs(estimatedUSCalibrationParameters[1] - trueUSCalibrationParameters[1])<TRANSLATION_EPS &&
            fabs(estimatedUSCalibrationParameters[2] - trueUSCalibrationParameters[2])<TRANSLATION_EPS &&
            anglesOK &&
            fabs(estimatedUSCalibrationParameters[6] - trueUSCalibrationParameters[6])<SCALE_EPS &&
            fabs(estimatedUSCalibrationParameters[7] - trueUSCalibrationParameters[7])<SCALE_EPS);
  }
}
