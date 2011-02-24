#include <fstream>
#include <time.h>
#include "SinglePointTargetUSCalibrationParametersEstimator.h"
#include "RANSAC.h"


bool loadTransformations(const std::string & transformationsFileName, 
                         std::vector<lsqrRecipes::Frame> &transformations);

bool loadImagePoints(const std::string &imagePointsFileName, 
                     size_t numberOfPoints, 
                     std::vector<lsqrRecipes::Point2D> &imagePoints);

bool loadSpatialPoints(const std::string &spatialPointsFileName, 
                       size_t numberOfPoints, 
                       std::vector<lsqrRecipes::Point3D> &spatialPoints);

int saveResults(const std::vector<double> &estimatedUSCalibrationParameters,
                const std::vector<lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator::DataType> &data,
                const std::vector<bool> &consensusSet,
                double percentageOfDataUsed,
                char *outputFileName);

std::string getCurrentDateTime(const char* format);

/**
 * Given the US-probe transformations, the corresponding 2D image coordinates and
 * 3D pointer coordinates compute the US calibration, write the result to the 
 * command line console and save in an xml file using the relevant IGSTK 
 * (www.igstk.org) tags for an affine transformation.
 *
 * @author Ziv Yaniv (zivy@isis.georgetown.edu)
 */
int main(int argc, char *argv[])
{
  typedef lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator::DataType DataType; 
  
  std::vector<lsqrRecipes::Frame> transformations;
  std::vector<lsqrRecipes::Point2D> imagePoints;
  std::vector<lsqrRecipes::Point3D> spatialPoints;
  DataType dataElement;
  std::vector< DataType > data;
  std::vector<double> estimatedUSCalibrationParameters;
  size_t i, n;

  if(argc != 5) {
    std::cerr<<"Unexpected number of command line parameters.\n";
    std::cerr<<"Usage: \n";
    std::cerr<<"\t"<<argv[0]<<" transformationsFileName 2DPointsFileName 3DPointsFileName outputFileName\n";
    return EXIT_FAILURE;
  }

  if(!loadTransformations(argv[1], transformations)) {
    std::cerr<<"Failed to load transformations file.\n";
    return EXIT_FAILURE;
  }
  if(!loadImagePoints(argv[2], transformations.size(), imagePoints)) {
    std::cerr<<"Failed to load image points file.\n";
    return EXIT_FAILURE;
  }
  if(!loadSpatialPoints(argv[3], transformations.size(), spatialPoints)) {
    std::cerr<<"Failed to load pointer points file.\n";
    return EXIT_FAILURE;
  }

  n = transformations.size();
  for(i=0; i<n; i++) {
    dataElement.T2 = transformations[i];
    dataElement.q = imagePoints[i];
    dataElement.p = spatialPoints[i];
    data.push_back(dataElement);
  }

  double maxDistanceBetweenPoints = 2.0;
  lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator 
    usCalibration(maxDistanceBetweenPoints);
  //usCalibration.setLeastSquaresType(lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator::ANALYTIC);

                          //estimate using the RANSAC algorithm
  double desiredProbabilityForNoOutliers = 0.999;
  double percentageOfDataUsed;
  std::vector<bool> consensusSet;
  percentageOfDataUsed = 
    lsqrRecipes::RANSAC< DataType, double>::compute(estimatedUSCalibrationParameters, 
                                                    &usCalibration, 
                                                    data, 
                                                    desiredProbabilityForNoOutliers,
                                                    &consensusSet);
  
  if(estimatedUSCalibrationParameters.size() == 0) {
   std::cout<<"FAILED CALIBRATION, possibly degenerate configuration\n\n\n";
   return EXIT_FAILURE;
  }
  else 
    return saveResults(estimatedUSCalibrationParameters, data, consensusSet, percentageOfDataUsed, argv[4]);
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
int saveResults(const std::vector<double> &estimatedUSCalibrationParameters,
                const std::vector<lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator::DataType> &data,
                const std::vector<bool> &consensusSet,
                double percentageOfDataUsed,
                char *outputFileName)
{          //should never happen, but just in case           
  if(estimatedUSCalibrationParameters.size() == 0)
   return EXIT_FAILURE;

  std::cout<<"Percentage of data used in estimate: "<<percentageOfDataUsed<<"\n";
  std::cout<<"t3[x,y,z]:\n";
  std::cout<<"\t["<<estimatedUSCalibrationParameters[0]<<", "<<estimatedUSCalibrationParameters[1];
  std::cout<<", "<<estimatedUSCalibrationParameters[2]<<"]\n";
  std::cout<<"omega[z,y,x]:\n";
  std::cout<<"\t["<<estimatedUSCalibrationParameters[3]<<", "<<estimatedUSCalibrationParameters[4];
  std::cout<<", "<<estimatedUSCalibrationParameters[5]<<"]\n";
  std::cout<<"m[x,y]:\n";
  std::cout<<"\t["<<estimatedUSCalibrationParameters[6];
  std::cout<<", "<<estimatedUSCalibrationParameters[7]<<"]\n";

  std::vector<double> errors;
  double minErr, maxErr, meanErr;

  std::vector<lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator::DataType> dataUsed;
  size_t i,n;
  n= data.size();
  for(i=0; i<n; i++) {
    if(consensusSet[i])
     dataUsed.push_back(data[i]);
  }

  lsqrRecipes::CalibratedPointerTargetUSCalibrationParametersEstimator::getDistanceStatistics(estimatedUSCalibrationParameters, dataUsed,
                                                                                               errors, minErr, maxErr, meanErr);
  double resNorm = 0.0;
  for(size_t z=0; z<errors.size(); z++)
  {
   resNorm+=(errors[z]*errors[z]);
  }
  std::cout<<"sum of squared errors: "<<resNorm<<"\n";
  std::cout<<"max, min, mean error: "<<maxErr<<", "<<minErr<<", "<<meanErr<<"\n\n\n";
                //write using the IGSTK (www.igstk.org) xml format tags
  std::ofstream out(outputFileName);
  if(!out.is_open())
   return EXIT_FAILURE;

            //ten digits after the decimal point should be enough to retain accuracy
            //in ASCII format
  const std::streamsize PRECISION = 10;
  out.precision(PRECISION);
  out.setf(std::ios::fixed,std::ios::floatfield); 

  out<<"<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n\n\n\n";
  out<<"<precomputed_transform>\n\n";
  out<<"\t<description>\n";
  out<<"\t"<<"US calibration - Pointer"<<"\n";  
  out<<"\t</description>\n\n";
  out<<"\t<computation_date>\n";
  out<<"\t"<<getCurrentDateTime( "%Y %b %d %H:%M:%S" )<<"\n";  
  out<<"\t</computation_date>\n\n";
  out<<"\t <transformation estimation_error=\"";
  out<<meanErr<<"\">\n";
  out<<"\t"<<estimatedUSCalibrationParameters[8]<<"\t"<<estimatedUSCalibrationParameters[11]<<"\t"<<estimatedUSCalibrationParameters[14]
     <<"\t"<< estimatedUSCalibrationParameters[0]<<"\n";
  out<<"\t"<<estimatedUSCalibrationParameters[9]<<"\t"<<estimatedUSCalibrationParameters[12]<<"\t"<<estimatedUSCalibrationParameters[15]
     <<"\t"<<estimatedUSCalibrationParameters[1]<<"\n";
  out<<"\t"<<estimatedUSCalibrationParameters[10]<<"\t"<<estimatedUSCalibrationParameters[13]<<"\t"<<estimatedUSCalibrationParameters[16]
     <<"\t"<< estimatedUSCalibrationParameters[2]<<"\n";
  out<<"\t</transformation>\n\n";
  out<<"</precomputed_transform>\n";
  out.close();
  return EXIT_SUCCESS;
}
/*****************************************************************************/
std::string getCurrentDateTime(const char* format)
{
  char buf[1024];
  time_t t;
  time(&t);
  strftime(buf, sizeof(buf), format, localtime(&t));
  return std::string(buf);
}
