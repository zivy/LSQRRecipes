#include "Frame.h"
#include "AbsoluteOrientationParametersEstimator.h"
#include "RandomNumberGenerator.h" 

bool reportAndCheckResults(const std::string &title,
                           const double distanceThreshold,
                           const std::vector<double> &estimatedTransformationParameters,
                           const std::vector< std::pair<lsqrRecipes::Point3D,lsqrRecipes::Point3D> > &evaluationData);

/*
 * Test the absolute orientation estimator's methods.
 */
int main(int argc, char *argv[])
{
  typedef std::pair<lsqrRecipes::Point3D,lsqrRecipes::Point3D> DataType;

	int i, pairNum = 10;
	double bounds = 100.0; //point coordinates are in [-100,100]
	double maxTranslation = 1000.0; //translations are in [-1000,1000]
	std::vector< DataType > estimationData, cleanEstimationData, targetData;
  double noiseSigma = 3.0; //noise is distributed IID ~N(0,noisSigma)
	DataType pointPair, outlierPair;
  std::vector<double> knownTransformationParameters, 
                      estimatedTransformationParameters;
  lsqrRecipes::Frame knownTransformation;

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

                      //create data without noise, data with noise and targets  
	for(i=0; i<pairNum; i++) {
              //data for estimation
		pointPair.first[0] = random.uniform(-bounds,bounds);
		pointPair.first[1] = random.uniform(-bounds,bounds);
		pointPair.first[2] = random.uniform(-bounds,bounds);
		knownTransformation.apply(pointPair.first,pointPair.second);
    cleanEstimationData.push_back(pointPair);
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
  outlierPair.first[0] = random.uniform(-bounds,bounds);
  outlierPair.first[1] = random.uniform(-bounds,bounds);
  outlierPair.first[2] = random.uniform(-bounds,bounds);
  knownTransformation.apply(outlierPair.first, outlierPair.second);
  outlierPair.second[0] +=10*noiseSigma;

	lsqrRecipes::AbsoluteOrientationParametersEstimator aopEstimator(1.0);
      //this threshold is used by the reportAndCheckResults() method which checks
      //that the maximal TRE for the random targets is less than this value
  double distanceThreshold = 3*noiseSigma;

	            //the known transformation
	std::cout<<"Known transformation:\n"<<knownTransformation;
           
		             //estimate the transformation parameters using only three point pairs
	aopEstimator.estimate(cleanEstimationData, estimatedTransformationParameters);
  if(!reportAndCheckResults("Estimation using three point pairs, no noise:",
                             distanceThreshold,
                             estimatedTransformationParameters,
                             targetData))
    return EXIT_FAILURE;
	      //least squares estimate
	aopEstimator.leastSquaresEstimate(estimationData, 
                                    estimatedTransformationParameters);
  if(!reportAndCheckResults("Least squares estimate, noisy data:",
                             distanceThreshold,
                             estimatedTransformationParameters,
                             targetData))
    return EXIT_FAILURE;
  
  if(!aopEstimator.agree(knownTransformationParameters, cleanEstimationData[0]) ||
      aopEstimator.agree(knownTransformationParameters, outlierPair))
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}

/*****************************************************************************/

bool reportAndCheckResults(const std::string &title,
                           const double distanceThreshold,
                           const std::vector<double> &estimatedTransformationParameters,
                           const std::vector< std::pair<lsqrRecipes::Point3D,lsqrRecipes::Point3D> > &evaluationData)
{
  std::cout<<"\n"<<title<<"\n\n";

  if(estimatedTransformationParameters.size() == 0) {
		std::cout<<"DEGENERATE CONFIGURATION\n\n";
    return false;
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
    std::cout<<"Estimated transformation:\n"<<estimatedTransformation<<"\n";
    unsigned int n = evaluationData.size();
    std::vector< std::pair<lsqrRecipes::Point3D,lsqrRecipes::Point3D> >::const_iterator it, end;
    end = evaluationData.end();
    double maxDistance = 0.0;
    for(it = evaluationData.begin(); it!=end; it++) {
      lsqrRecipes::Point3D pTransformed;
      estimatedTransformation.apply((*it).first,pTransformed);
      double distance = (pTransformed - (*it).second).l2Norm();
      if(distance>maxDistance)
        maxDistance = distance;
    }
    std::cout<<"Maximal TRE to random targets is: "<<maxDistance<<"\n";
    if(maxDistance<distanceThreshold)
      return true;
    else
      return false;
  }
        //get rid of the warnning about not all control paths returning a value
  return false;
}
