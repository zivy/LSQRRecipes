#include "Frame.h"
#include "AbsoluteOrientationParametersEstimator.h"
#include "RandomNumberGenerator.h" 
#include "RANSAC.h"

void reportResults(const std::string &title,
                   const std::vector<double> &estimatedTransformationParameters,
                   const std::vector< std::pair<lsqrRecipes::Point3D,lsqrRecipes::Point3D> > &estimationData,
                   const std::vector< std::pair<lsqrRecipes::Point3D,lsqrRecipes::Point3D> > &targetData);

/*
 * Example of using the exhaustive search method that is part of the RANSAC
 * implementation. For absolute orientation with a small number of points this
 * is computationally reasonable. The program estimates the transformation using
 * a least squares method and a RANSAC wrapper for the same least squares method.
 * The program then prints the estimated/known transformations and associated
 * FREs and TREs. Look at the FRE to identify the outlying point(s).
 */
int main(int argc, char *argv[])
{
  typedef std::pair<lsqrRecipes::Point3D,lsqrRecipes::Point3D> DataType;

  unsigned int inliers = 4;
  unsigned int outliers = 1;
	unsigned int i;
	double bounds = 100.0; //point coordinates are in [-100,100]
	double maxTranslation = 1000.0; //translations are in [-1000,1000]
	std::vector< DataType > estimationData, targetData;
  double noiseSigma = 1.0; //noise is distributed IID ~N(0,noisSigma)
	DataType pointPair;
  std::vector<double> knownTransformationParameters, 
                      estimatedTransformationParameters;
  lsqrRecipes::Frame knownTransformation, estimatedTransformation;

	lsqrRecipes::RandomNumberGenerator random;

	                          //create a random transformation
	            //random unit quaternion
  double qx,qy,qz;
	qx = random.uniform(0.0,1.0);
	qy = random.uniform(0.0,sqrt(1.0-qx*qx));
	qz = random.uniform(0.0,sqrt(1.0-qx*qx-qy*qy));
  knownTransformationParameters.push_back(sqrt(1.0 - qx*qx - qy*qy - qz*qz));
  knownTransformationParameters.push_back(qx);
  knownTransformationParameters.push_back(qy);
  knownTransformationParameters.push_back(qz);
               //random translation
  knownTransformationParameters.push_back(random.uniform(-maxTranslation,maxTranslation));
  knownTransformationParameters.push_back(random.uniform(-maxTranslation,maxTranslation));
  knownTransformationParameters.push_back(random.uniform(-maxTranslation,maxTranslation));

  knownTransformation.setRotationQuaternion(knownTransformationParameters[0],
                                            knownTransformationParameters[1], 
                                            knownTransformationParameters[2],
                                            knownTransformationParameters[3]);
  knownTransformation.setTranslation(knownTransformationParameters[4],
                                     knownTransformationParameters[5], 
                                     knownTransformationParameters[6]);

                      //create data and targets  
	for(i=0; i<inliers; i++) {
              //inliers
		pointPair.first[0] = random.uniform(-bounds,bounds);
		pointPair.first[1] = random.uniform(-bounds,bounds);
		pointPair.first[2] = random.uniform(-bounds,bounds);
		knownTransformation.apply(pointPair.first,pointPair.second);
    pointPair.second[0] += random.normal(noiseSigma);
		pointPair.second[1] += random.normal(noiseSigma);
		pointPair.second[2] += random.normal(noiseSigma);
		estimationData.push_back(pointPair);
             //data for evaluation, our random targets
		pointPair.first[0] = random.uniform(-bounds,bounds);
		pointPair.first[1] = random.uniform(-bounds,bounds);
		pointPair.first[2] = random.uniform(-bounds,bounds);
		knownTransformation.apply(pointPair.first,pointPair.second);
    targetData.push_back(pointPair);
	}
             //add our outliers
	for(i=0; i<outliers; i++) {              
		pointPair.first[0] = random.uniform(-bounds,bounds);
		pointPair.first[1] = random.uniform(-bounds,bounds);
		pointPair.first[2] = random.uniform(-bounds,bounds);
		knownTransformation.apply(pointPair.first,pointPair.second);
    pointPair.second[0] += random.normal(noiseSigma);
		pointPair.second[1] += random.normal(noiseSigma);
		pointPair.second[2] += 5.0;
		estimationData.push_back(pointPair);
  }

  reportResults("known transformation:", knownTransformationParameters,
                estimationData, targetData);


       //a pair of points is consistent with a transformation if
       // ||T*first - second|| < distanceThreshold
  double distanceThreshold = 2*noiseSigma;
	lsqrRecipes::AbsoluteOrientationParametersEstimator aopEstimator(distanceThreshold);

                 //standard least squares estimate
  aopEstimator.leastSquaresEstimate(estimationData, 
                                    estimatedTransformationParameters);
  reportResults("least squares estimate:", estimatedTransformationParameters,
                estimationData, targetData);

                 //exhaustive search estimate
  std::vector<bool>consensusSet;
  double percentageOfDataUsed;
  percentageOfDataUsed = 
    lsqrRecipes::RANSAC< std::pair<lsqrRecipes::Point3D,lsqrRecipes::Point3D>, double>::compute(estimatedTransformationParameters, 
                                                      &aopEstimator, 
                                                      estimationData, &consensusSet);
  reportResults("exhaustive search estimate:", estimatedTransformationParameters,
                estimationData, targetData);
  std::cout<<"Fiducials used in final estimate: ";
  for(std::vector<bool>::const_iterator it = consensusSet.begin(); it!=consensusSet.end(); it++)
    std::cout<<(*it)<<" ";
  std::cout<<"\n";

  return EXIT_SUCCESS;
}

/*****************************************************************************/
void reportResults(const std::string &title,
                   const std::vector<double> &estimatedTransformationParameters,
                   const std::vector< std::pair<lsqrRecipes::Point3D,lsqrRecipes::Point3D> > &estimationData,
                   const std::vector< std::pair<lsqrRecipes::Point3D,lsqrRecipes::Point3D> > &targetData)
{
  std::cout<<"\n"<<title<<"\n\n";

  if(estimatedTransformationParameters.size() == 0) {
		std::cout<<"DEGENERATE CONFIGURATION\n\n";
  }
	else {
    lsqrRecipes::Frame estimatedTransformation;
    estimatedTransformation.setRotationQuaternion(estimatedTransformationParameters[0],
                                                  estimatedTransformationParameters[1], 
                                                  estimatedTransformationParameters[2],
                                                  estimatedTransformationParameters[3]);
    estimatedTransformation.setTranslation(estimatedTransformationParameters[4],
                                           estimatedTransformationParameters[5], 
                                           estimatedTransformationParameters[6]);
    std::cout<<"Transformation:\n"<<estimatedTransformation<<"\n";
             //look at the FRE
    std::vector< std::pair<lsqrRecipes::Point3D,lsqrRecipes::Point3D> >::const_iterator it, end;
    end = estimationData.end();
    std::cout<<"Fiducial registration errors: ";
    for(it = estimationData.begin(); it!=end; it++) {
      lsqrRecipes::Point3D pTransformed;
      estimatedTransformation.apply((*it).first,pTransformed);
      double distance = (pTransformed - (*it).second).l2Norm();
      std::cout<<distance<<" ";
    }
    std::cout<<"\n";
             //look at the TRE
    end = targetData.end();
    std::cout<<"Target registration errors: ";
    for(it = targetData.begin(); it!=end; it++) {
      lsqrRecipes::Point3D pTransformed;
      estimatedTransformation.apply((*it).first,pTransformed);
      double distance = (pTransformed - (*it).second).l2Norm();
      std::cout<<distance<<" ";
    }
    std::cout<<"\n";
  }
}
