#include <fstream>
#include "Ray3D.h"
#include "RandomNumberGenerator.h"
#include "RayIntersectionParametersEstimator.h"
#include "RANSAC.h"

/**
 * Generate rays that intersect at a single point, additive zero mean Gaussian 
 * noise, and outliers.
 * @param numInliers How many rays are inliers.
 * @param numOutliers How many rays are outliers.
 * @param outlierDistance Threshold defining outliers, the distance from the ray
 *                        to the intersection point is larger than this.
 * @param data The rays are added to the end of this vector.
 * @param intersectionParameters [a], 3D point. 
 */
void generateData(unsigned int numInliers, unsigned int numOutliers,
                  double outlierDistance,
                  std::vector< lsqrRecipes::Ray3D > &data,
                  std::vector<double> &intersectionParameters);

/**
 * Save an open inventor ASCII scene file showing the rays identified as
 * inliers/outliers and the estimated intersection point. The format is 
 * compatible with the viewers available in 
 * the open source coin3D toolkit (www.coin3d.org).
 *
 * @param outputFileName Write the scene to this file.
 * @param data The 3D points.
 * @param rayLength The length of the ray segment to draw.
 * @param estimatedPlaneParameters [n,a], plane normal and point on plane. Plane 
 *                                 is the set of points p such that n^T(p-a)=0.
 * @param parameterEstimator The plane parameter estimator whose agree() method
 *                           is used to identify inliers.
 */
void saveOIVFile(std::string &outputFileName,
                 std::vector< lsqrRecipes::Ray3D > &data,
                 double rayLength,
                 std::vector<double> &estimatedIntersectionParameters,
                 lsqrRecipes::RayIntersectionParametersEstimator 
                 &parameterEstimator);

/**
 * Given the number of outliers and inliers generate a random dataset 
 * accordingly. Then estimate the intersection point using a least squares 
 * estimate and the RANSAC algorithm. Compare the results
 * to the known intersection point. Two open inventor scene files are written, 
 * showing the least squares and RANSAC estimates. The rays' origin is marked
 * by a sphere. Rays that agree with the estimated model are green, otherwise 
 * they are red.
 *
 * @author Ziv Yaniv 
 */
int main(int argc, char *argv[])
{
  const unsigned int INLIERS = 9;
  const unsigned int OUTLIERS = 1;
        //length of ray we use for drawing
  double l, rayLength;
  std::string leastSquaresOutputFileName = "leastSquaresRayIntersection.iv";
  std::string ransacOutputFileName = "RANSACRayIntersection.iv";

  std::vector<lsqrRecipes::Ray3D> data;
  std::vector<double> trueIntersectionParameters, intersectionParameters; 
                          //double the noise level sigma in GenerateData()
  double maximalDistanceFromRay = 4.0;   
  
  generateData(INLIERS, OUTLIERS, maximalDistanceFromRay, data, 
               trueIntersectionParameters);
          //compute ray length we use for drawing
  rayLength = 0.0;
  lsqrRecipes::Point3D truePoint;
  truePoint[0] = trueIntersectionParameters[0];
  truePoint[1] = trueIntersectionParameters[1];
  truePoint[2] = trueIntersectionParameters[2];
  for(size_t i=0; i<data.size(); i++) {
    if((l=(truePoint - data[i].p).l2Norm())>rayLength)
      rayLength = l;
  }
  rayLength*=2;

  std::cout<<"Known intersection parameters [a]\n\t [ ";
  std::cout<<trueIntersectionParameters[0]<<", "<<trueIntersectionParameters[1];
  std::cout<<", "<<trueIntersectionParameters[2]<<"]\n\n";

                        //create and initialize the parameter estimator
  
  lsqrRecipes::RayIntersectionParametersEstimator ripEstimator(maximalDistanceFromRay);
  ripEstimator.leastSquaresEstimate(data, intersectionParameters);
  if(intersectionParameters.empty())
    std::cout<<"Least squares estimate failed, degenerate configuration?\n";
  else
  {
	  std::cout<<"Least squares intersection parameters: [a]\n\t [ ";
    std::cout<<intersectionParameters[0]<<", "<<intersectionParameters[1];
    std::cout<<", "<<intersectionParameters[2]<<"]\n\n";
             //save scene file
    saveOIVFile(leastSquaresOutputFileName, data, rayLength, 
                intersectionParameters, ripEstimator);
  }
                          //estimate using the RANSAC algorithm
  double desiredProbabilityForNoOutliers = 0.999;
  double percentageOfDataUsed;
  percentageOfDataUsed = 
    lsqrRecipes::RANSAC<lsqrRecipes::Ray3D, double>::compute(intersectionParameters, 
                                                              &ripEstimator, 
                                                              data, 
                                                              desiredProbabilityForNoOutliers);
  if(intersectionParameters.empty())
    std::cout<<"RANSAC estimate failed, degenerate configuration?\n";
  else
  {
	  std::cout<<"RANSAC intersection parameters: [a]\n\t [ ";
    std::cout<<intersectionParameters[0]<<", "<<intersectionParameters[1];
    std::cout<<", "<<intersectionParameters[2]<<"]\n\n";
             //save scene file
    saveOIVFile(ransacOutputFileName, data, rayLength, 
                intersectionParameters, ripEstimator);
  }
  return EXIT_SUCCESS;
}

/*****************************************************************************/

void generateData(unsigned int numInliers, unsigned int numOutliers,
                  double outlierDistance,
                  std::vector< lsqrRecipes::Ray3D > &data,
                  std::vector<double> &intersectionParameters)
{
  unsigned int i;
  lsqrRecipes::Ray3D ray;
  lsqrRecipes::RayIntersectionParametersEstimator rip(outlierDistance);

         //noise is distributed IID ~N(0,noiseSigma)
  const double NOISE_SIGMA = 2.0; 
         //all data is generated in [-maxRange,maxRange]X[-maxRange,maxRange]X[-maxRange,maxRange]
  double maxRange =1000.0;

  lsqrRecipes::Point3D knownIntersectionPoint;  

  lsqrRecipes::RandomNumberGenerator random;
             //randomly select intersection point
  knownIntersectionPoint[0] = random.uniform(-maxRange,maxRange);
  knownIntersectionPoint[1] = random.uniform(-maxRange,maxRange);
  knownIntersectionPoint[2] = random.uniform(-maxRange,maxRange);
  intersectionParameters.clear();
  intersectionParameters.push_back(knownIntersectionPoint[0]);
  intersectionParameters.push_back(knownIntersectionPoint[1]);
  intersectionParameters.push_back(knownIntersectionPoint[2]);

             //randomly create rays approximately going through the given point
  for(i=0; i<numInliers; i++) {
    ray.p[0] = random.uniform(-maxRange,maxRange);
    ray.p[1] = random.uniform(-maxRange,maxRange);
    ray.p[2] = random.uniform(-maxRange,maxRange);
    ray.n[0] = knownIntersectionPoint[0] + random.normal(NOISE_SIGMA) - ray.p[0];
    ray.n[1] = knownIntersectionPoint[1] + random.normal(NOISE_SIGMA) - ray.p[1];
    ray.n[2] = knownIntersectionPoint[2] + random.normal(NOISE_SIGMA) - ray.p[2];
    ray.n.normalize();
    data.push_back(ray);
  }
              //add outliers, via rejection method
  for(i=0; i<numOutliers; i++) {
    ray.p[0] = random.uniform(-maxRange,maxRange);
    ray.p[1] = random.uniform(-maxRange,maxRange);
    ray.p[2] = random.uniform(-maxRange,maxRange);
    ray.n[0] = random.uniform(-maxRange,maxRange) - ray.p[0];
    ray.n[1] = random.uniform(-maxRange,maxRange) - ray.p[1];
    ray.n[2] = random.uniform(-maxRange,maxRange) - ray.p[2];
    ray.n.normalize();
    if(rip.agree(intersectionParameters, ray))
     i--;
    else
     data.push_back(ray);
  }
}

/*****************************************************************************/

void saveOIVFile(std::string &outputFileName,
                 std::vector< lsqrRecipes::Ray3D > &data,
                 double rayLength,
                 std::vector<double> &estimatedIntersectionParameters,
                 lsqrRecipes::RayIntersectionParametersEstimator 
                 &parameterEstimator)
{
  double sphereRadius = 50;
               //outliers are metallic red 
  std::string outlierMaterial = "Material {\n";
  outlierMaterial+="\tambientColor 1.0 0.0 0.0\n";
  outlierMaterial+="\tdiffuseColor 0.27 0.15 0.0\n";
  outlierMaterial+="\tspecularColor 1.0 0.0 0.0\n";
  outlierMaterial+="}\n";

               //inliers are metallic green 
  std::string inlierMaterial = "\tMaterial {\n";
  inlierMaterial+="\t\tambientColor 0.0 1.0 0.0\n";
  inlierMaterial+="\t\tdiffuseColor 0.0 0.27 0.15\n";
  inlierMaterial+="\t\tspecularColor 0.0 1.0 0.0\n";
  inlierMaterial+="\t}\n";

             //intersection point is gray
  std::string intersectionMaterial = "\tMaterial {\n";
  intersectionMaterial+="\t\tambientColor 0.5 0.5 0.5\n";
  intersectionMaterial+="\t\tdiffuseColor 0.2 0.2 0.2\n";
  intersectionMaterial+="\t\tspecularColor 0.8 0.8 0.8\n";
  intersectionMaterial+="\t}\n";

  std::ofstream out( outputFileName.c_str() );
  if(!out.is_open())
   return;

  out<<"#Inventor V2.1 ascii\n\n";
               //go over all the rays and write the outliers and inliers to the 
               //scene graph
  for(unsigned int i=0; i<data.size(); i++) {
    out<<"Separator {\n";
    if(parameterEstimator.agree(estimatedIntersectionParameters, data[i]))
      out<<inlierMaterial;    
    else
      out<<outlierMaterial;
             //write the sphere at the ray's origin
    out<<"\tSeparator {\n";
    out<<"\tTransform {\n";
    out<<"\t\ttranslation "<<(data[i]).p[0]<<" "<<(data[i]).p[1]<<" "<<(data[i]).p[2]<<"\n";
    out<<"\t}\n";
    out<<"\tSphere {\n";
    out<<"\t\tradius  "<<sphereRadius<<"\n";
    out<<"\t}\n";
    out<<"}\n";
    out<<"\tCoordinate3 {\n";
    out<<"\t\tpoint ["<<(data[i]).p[0]<<" "<<(data[i]).p[1]<<" "<<(data[i]).p[2]<<", ";
    out<<(data[i]).p[0]+(data[i]).n[0]*rayLength<<" ";
    out<<(data[i]).p[1]+(data[i]).n[1]*rayLength<<" ";
    out<<(data[i]).p[2]+(data[i]).n[2]*rayLength<<"]";
    out<<"\t\t}\n";
    out<<"\t\tLineSet { numVertices [-1] }\n";
    out<<"\t}\n";
  }

  out<<"\tSeparator {\n";
  out<<intersectionMaterial;
  out<<"\tTransform {\n";
  out<<"\t\ttranslation "<<estimatedIntersectionParameters[0]<<" ";
  out<<estimatedIntersectionParameters[1]<<" ";
  out<<estimatedIntersectionParameters[2]<<"\n";
  out<<"\t}\n";
  out<<"\tSphere {\n";
  out<<"\t\tradius  "<<sphereRadius<<"\n";
  out<<"\t}\n";
  out<<"}\n";
  out.close();
}
