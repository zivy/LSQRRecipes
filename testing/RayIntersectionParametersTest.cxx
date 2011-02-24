#include <stdlib.h>
#include "RandomNumberGenerator.h"
#include "Ray3D.h"
#include "RayIntersectionParametersEstimator.h"


/*
 * Test the plane estimator's methods. We arbitrarily choose to work with 3D
 * planes even though the tested code works for any dimension.
 */
int main( int argc, char *argv[] )
{
  std::vector<lsqrRecipes::Ray3D> rayData;
  std::vector<lsqrRecipes::Ray3D> noNoise;

  lsqrRecipes::Ray3D ray;
  unsigned int i;
          //number of rays that intersect at the same point
	const unsigned int NUM_RAYS = 10; 
         //noise is distributed IID ~N(0,noiseSigma)
  const double NOISE_SIGMA = 20.0; 
         //all data is generated in [-maxRange,maxRange]X[-maxRange,maxRange]X[-maxRange,maxRange]
  double maxRange =1000.0;

  lsqrRecipes::Point3D knownIntersectionPoint;  

  lsqrRecipes::RandomNumberGenerator random;
             //randomly select intersection point
  knownIntersectionPoint[0] = random.uniform(-maxRange,maxRange);
  knownIntersectionPoint[1] = random.uniform(-maxRange,maxRange);
  knownIntersectionPoint[2] = random.uniform(-maxRange,maxRange);

             //randomly create rays intersecting at approximatly the given point
  for(i=0; i<NUM_RAYS; i++) {
    ray.p[0] = random.uniform(-maxRange,maxRange);
    ray.p[1] = random.uniform(-maxRange,maxRange);
    ray.p[2] = random.uniform(-maxRange,maxRange);
    ray.n[0] = knownIntersectionPoint[0] + random.normal(NOISE_SIGMA) - ray.p[0];
    ray.n[1] = knownIntersectionPoint[1] + random.normal(NOISE_SIGMA) - ray.p[1];
    ray.n[2] = knownIntersectionPoint[2] + random.normal(NOISE_SIGMA) - ray.p[2];
    ray.n.normalize();
    rayData.push_back(ray);
  }
            //two rays intersecting at the given point   
  for(i=0; i<2; i++) {
    ray.p[0] = random.uniform(-maxRange,maxRange);
    ray.p[1] = random.uniform(-maxRange,maxRange);
    ray.p[2] = random.uniform(-maxRange,maxRange);
    ray.n[0] = knownIntersectionPoint[0] - ray.p[0];
    ray.n[1] = knownIntersectionPoint[1] - ray.p[1];
    ray.n[2] = knownIntersectionPoint[2] - ray.p[2];
    ray.n.normalize();
    noNoise.push_back(ray);
  }
  
  

          //2. Test the ray intersection code
          //(a) test the agree() method
          //(b) intersection of two rays, no noise
	        //(c) least squares estimate with noisy data

	std::vector<double> estimatedIntersectionPoint;
  double maxDistanceToRay = 0.5;
	lsqrRecipes::RayIntersectionParametersEstimator ripEstimator(maxDistanceToRay);
	
	            //The known intersection point
  std::cout<<"Known intersection point [x,y,z]:\n\t [ "<<knownIntersectionPoint[0]<<", ";
  std::cout<<knownIntersectionPoint[1]<<", "<<knownIntersectionPoint[2]<<" ]\n\n";
             //test the agree() method
  estimatedIntersectionPoint.push_back(knownIntersectionPoint[0]);
  estimatedIntersectionPoint.push_back(knownIntersectionPoint[1]);
  estimatedIntersectionPoint.push_back(knownIntersectionPoint[2]);
  if(!ripEstimator.agree(estimatedIntersectionPoint, noNoise[0]))
    return EXIT_FAILURE;

             //Compute the intersection point using two rays, no noise
  ripEstimator.estimate(noNoise, estimatedIntersectionPoint);
	if(estimatedIntersectionPoint.size() == 0)
		std::cout<<"Intersection point using two rays [x,y,z]: DEGENERATE CONFIGURATION\n\n";
	else {
		std::cout<<"Intersection point using two rays [x,y,z], no noise:\n\t [ ";
    std::cout<<estimatedIntersectionPoint[0]<<", "<<estimatedIntersectionPoint[1]<<", ";
		std::cout<<estimatedIntersectionPoint[2]<<" ]\n\n";
         //check that the known and estimated points are the same, they are 
         //expected to be the same as there is no noise
    lsqrRecipes::Point3D tmp;
    tmp[0] = estimatedIntersectionPoint[0];
    tmp[1] = estimatedIntersectionPoint[1];
    tmp[2] = estimatedIntersectionPoint[2];
             //use same distance threshold as used for determining if a point is
             //on the ray or not    
    if((tmp-knownIntersectionPoint).l2Norm()> maxDistanceToRay)
     return EXIT_FAILURE;
	}
             //least squares estimate of the intersection point
  ripEstimator.leastSquaresEstimate(rayData, estimatedIntersectionPoint);
	if(estimatedIntersectionPoint.size() == 0)
		std::cout<<"Intersection point using least squares [x,y,z]: DEGENERATE CONFIGURATION\n\n";
	else {
   std::cout<<"Intersection point using least squares [x,y,z]:\n\t [ ";
   std::cout<<estimatedIntersectionPoint[0]<<", "<<estimatedIntersectionPoint[1]<<", ";
	 std::cout<<estimatedIntersectionPoint[2]<<" ]\n\n";
	}
  return EXIT_SUCCESS;
}
