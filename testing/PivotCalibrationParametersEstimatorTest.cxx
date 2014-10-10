#include <fstream>
#include <stdlib.h>
#include "PivotCalibrationParametersEstimator.h"

using namespace lsqrRecipes;

int main(int argc, char *argv[])
{
  if(argc != 2) {
    std::cerr<<"Unexpected number of command line parameters.\n";
    std::cerr<<"Usage: \n";
    std::cerr<<"\t"<<argv[0]<<" pivotCalibrationDataFileName\n";
    return EXIT_FAILURE;
  }

  double x, y, z, qx, qy, qz, qs;
  std::vector< Frame > transformations;
  Frame f;
  unsigned int i;
  const unsigned int PARAMETER_NUM = 6;

               //load the data, each row in the file is expected to have the
               //following format x y z qx qy qz qs
  std::ifstream in;

  in.open(argv[1]);
  if(!in.is_open())
    return EXIT_FAILURE;
  while(in>> x >> y >> z >> qx >> qy >> qz >> qs) {
    f.setRotationQuaternion(qs, qx, qy, qz);
    f.setTranslation(x, y, z);   
    transformations.push_back(f);
  }
  in.close();

  if(transformations.empty()) {
    std::cerr<<"Failed to load data from pivot calibration file.\n";
    return EXIT_FAILURE;
  }

      //configure the estimator
  double maxError = 1.0; //distance between the two points is at most 1mm
  std::vector<double> estimatedTranslations;
  PivotCalibrationEstimator pivotCalibration(maxError);

            //test the exact estimate method
  double knownExactTranslations[] = {-18.586, 1.98134, -157.439, 
                                      146.965, -62.0497, -1042.87};

  std::vector<Frame> minTransformations(3);
  minTransformations[0] = transformations[0];
  minTransformations[1] = transformations[static_cast<unsigned int>(transformations.size()/2.0)];
  minTransformations[2] = transformations[transformations.size()-1];
  
  pivotCalibration.estimate(minTransformations, estimatedTranslations);
  if(estimatedTranslations.empty())
    return EXIT_FAILURE;

  std::cout<<"Exact estimate\n";
  std::cout<<"--------------\n";
  std::cout<<"Known translations [DRF^t, W^t]:\n\t [ ";
  for(i=0; i<PARAMETER_NUM-1; i++)
    std::cout<<knownExactTranslations[i]<<", ";
  std::cout<<knownExactTranslations[PARAMETER_NUM-1]<<"]\n\n";
  std::cout<<"Estimated translations [DRF^t, W^t]:\n\t [ ";
  for(i=0; i<PARAMETER_NUM-1; i++)
    std::cout<<estimatedTranslations[i]<<", ";
  std::cout<<estimatedTranslations[PARAMETER_NUM-1]<<"]\n\n";

  bool succeedExact = true;
  for(i=0; i<PARAMETER_NUM; i++ )
    succeedExact = succeedExact && (fabs(estimatedTranslations[i] - knownExactTranslations[i])<maxError);
 
          //test the agree method
  for(i=0;i<minTransformations.size(); i++) {
    if(!pivotCalibration.agree(estimatedTranslations, minTransformations[i]))
      return EXIT_FAILURE;
  }


            //test the least squares estimation method
  double knownLSTranslations[] = {-17.7799, 1.1113, -156.865, 
                                   146.901, -62.9689, -1042.14};

  pivotCalibration.leastSquaresEstimate(transformations, estimatedTranslations);
  if(estimatedTranslations.empty())
    return EXIT_FAILURE;

  std::cout<<"Least Squares estimate\n";
  std::cout<<"----------------------\n";
  std::cout<<"Known translations [DRF^t, W^t]:\n\t [ ";
  for(i=0; i<PARAMETER_NUM-1; i++)
    std::cout<<knownLSTranslations[i]<<", ";
  std::cout<<knownLSTranslations[PARAMETER_NUM-1]<<"]\n\n";
  std::cout<<"Estimated translations [DRF^t, W^t]:\n\t [ ";
  for(i=0; i<PARAMETER_NUM-1; i++)
    std::cout<<estimatedTranslations[i]<<", ";
  std::cout<<estimatedTranslations[PARAMETER_NUM-1]<<"]\n\n";

  bool succeedLeastSquares = true;
  for(i=0; i<PARAMETER_NUM; i++ )
    succeedLeastSquares = succeedLeastSquares && (fabs(estimatedTranslations[i] - knownLSTranslations[i])<maxError);

  if(succeedExact && succeedLeastSquares)
    return EXIT_SUCCESS;
  return EXIT_FAILURE;
}

