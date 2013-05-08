#include <fstream>
#include "RANSAC.h"
#include "RandomNumberGenerator.h"
#include "SphereParametersEstimator.h"


/**
 * Generate points on a hypersphere with additive zero mean Gaussian noise and
 * outliers.
 * @param numInliers How many points are inliers.
 * @param numOutliers How many points are outliers.
 * @param outlierDistance Threshold defining outliers, points that are farther
 *                        than this distance from the sphere.
 * @param data The points are added to the end of this vector.
 * @param sphereParameters [c,r], sphere center and radius. Sphere is 
 *                        defined as the set of points p such that 
 *                        (p-c)^T(p-c)=r^2. 
 */
template<unsigned int dimension>
void generateData( unsigned int numInliers, unsigned int numOutliers, 
                   double outlierDistance,
                   std::vector< lsqrRecipes::Point<double,dimension> > &data,
                   std::vector<double> &sphereParameters );

/**
 * This function only works in 3D.
 *
 * Save an open inventor ASCII scene file showing the points identified as
 * inliers/outliers and the estimated sphere. The output file is 
 * "sphereEstimation.iv". The format is compatible with the viewers available in 
 * the open source coin3D toolkit (www.coin3d.org).
 *
 * @param outputFileName Write the scene to this file.
 * @param data The 3D points.
 * @param estimatedSphereParameters [c,r], sphere center and radius. Sphere is 
 *                                  defined as the set of points p such that 
 *                                 (p-c)^T(p-c)=r^2. 
 * @param parameterEstimator The sphere parameter estimator whose agree() method
 *                           is used to identify inliers.
 */
template<unsigned int dimension>
void saveOIVFile(std::string &outputFileName,
                 std::vector< lsqrRecipes::Point<double,dimension> > &data, 
                 std::vector<double> &estimatedSphereParameters,
                 lsqrRecipes::SphereParametersEstimator<dimension> 
                   & parameterEstimator);

/**
 * Given the hard coded dimension, and number of outliers and inliers generate
 * a random dataset accordingly. Then estimate the hypersphere parameter values
 * using a least squares estimate and the RANSAC algorithm. Compare the results
 * to the known hypersphere. Code is written for nD data except the 
 * visualization which is limited to 3D. If DIMENSION is set to three, then two
 * open inventor scene files are written, showing the least squares and RANSAC
 * estimates. Data points are colored spheres, those that agree with the 
 * estimated model are green, otherwise they are red.
 *
 * @author Ziv Yaniv 
 */
int main(int argc, char *argv[])
{
  const unsigned int DIMENSION = 3;
  const unsigned int INLIERS = 90;
  const unsigned int OUTLIERS = 10;
  std::string leastSquaresOutputFileName = "leastSquaresSphereEstimation.iv";
  std::string ransacOutputFileName = "RANSACSphereEstimation.iv";
  
  std::vector< lsqrRecipes::Point<double,DIMENSION> > data;
  std::vector<double> trueSphereParameters, sphereParameters; 
  double outlierDistance = 20.0;
  unsigned int i;
  lsqrRecipes::Vector<double, DIMENSION> tmp; 

  generateData<DIMENSION>( INLIERS, OUTLIERS, outlierDistance,
                           data, trueSphereParameters );

  std::cout<<"Known (hyper)sphere parameters [c,r]\n\t [ ";
  for( i=0; i<DIMENSION; i++ )
    std::cout<<trueSphereParameters[i]<<", ";
  std::cout<<trueSphereParameters[i]<<"]\n\n";
                        //create and initialize the parameter estimator
  double maximalDistanceFromSphere = 0.5;
  lsqrRecipes::SphereParametersEstimator<DIMENSION> sphereEstimator(maximalDistanceFromSphere,
                                                                    lsqrRecipes::SphereParametersEstimator<DIMENSION>::GEOMETRIC);
  sphereEstimator.leastSquaresEstimate( data, sphereParameters );
  if( sphereParameters.empty() )
    std::cout<<"Least squares estimate failed, degenerate configuration?\n";
  else {
	  std::cout<<"Least squares (hyper)sphere parameters: [c,r]\n\t [ ";
    for( i=0; i<DIMENSION; i++ )
      std::cout<<sphereParameters[i]<<", ";
    std::cout<<sphereParameters[i]<<"]\n\n";
              //distance between known hypersphere center and estimated one,
              //and difference between the two radii
    for( i=0; i<DIMENSION; i++ )
      tmp[i] = sphereParameters[i] - trueSphereParameters[i];
    std::cout<<"\t Distance between estimated and known sphere centers [0=correct]: ";
    std::cout<<tmp.l2Norm()<<"\n";
    std::cout<<"\t Difference between estimated and known sphere radius [0=correct]: ";
    std::cout<<fabs( sphereParameters[DIMENSION] - trueSphereParameters[DIMENSION] )<<"\n";
              //save scene file (works only in 3D)
    saveOIVFile<DIMENSION>( leastSquaresOutputFileName, data, sphereParameters, sphereEstimator );
  }
                          //estimate using the RANSAC algorithm
  double desiredProbabilityForNoOutliers = 0.999;
  double percentageOfDataUsed;
  percentageOfDataUsed = 
    lsqrRecipes::RANSAC< lsqrRecipes::Point<double, DIMENSION>, double>::compute(sphereParameters, 
                                                      &sphereEstimator, 
                                                      data, 
                                                      desiredProbabilityForNoOutliers);
  if( sphereParameters.empty() )
    std::cout<<"RANSAC estimate failed, degenerate configuration?\n";
  else
  {
	  std::cout<<"RANSAC (hyper)sphere parameters: [c,r]\n\t [ ";
    for( i=0; i<DIMENSION; i++ )
      std::cout<<sphereParameters[i]<<", ";
    std::cout<<sphereParameters[i]<<"]\n\n";
              //distance between known hypersphere center and estimated one,
              //and difference between the two radii
    for( i=0; i<DIMENSION; i++ )
      tmp[i] = sphereParameters[i] - trueSphereParameters[i];
    std::cout<<"\t Distance between estimated and known sphere centers [0=correct]: ";
    std::cout<<tmp.l2Norm()<<"\n";
    std::cout<<"\t Difference between estimated and known sphere radius [0=correct]: ";
    std::cout<<fabs( sphereParameters[DIMENSION] - trueSphereParameters[DIMENSION] )<<"\n";
            //save scene file (works only in 3D)
    saveOIVFile<DIMENSION>( ransacOutputFileName, data, sphereParameters, sphereEstimator );
  }
  return EXIT_SUCCESS;
}

/*****************************************************************************/

template<unsigned int dimension>
void generateData( unsigned int numInliers, unsigned int numOutliers, 
                   double outlierDistance,
                   std::vector< lsqrRecipes::Point<double,dimension> > &data,
                   std::vector<double> &sphereParameters )
{
  lsqrRecipes::Vector<double, dimension> tmp, noise;
  lsqrRecipes::Point<double, dimension> sphereCenter, randomPoint;
  double sphereRadius;
  double noiseStandardDeviation = 0.4; //noise standard deviation
  double coordinateMax = 1000.0;
  unsigned int i, j;

  lsqrRecipes::RandomNumberGenerator random;

  sphereParameters.clear();
         //generate points on random hypersphere
  for( i=0; i<dimension; i++ ) {
    sphereCenter[i] = random.uniform( -coordinateMax, coordinateMax );
  }
  sphereRadius = random.uniform(0.0, coordinateMax);

  for( i=0; i<dimension; i++ ) 
    sphereParameters.push_back( sphereCenter[i] );
  sphereParameters.push_back( sphereRadius );

               //generate inliers
  for( i=0; i<numInliers; i++ ) {
    for( j=0; j<dimension; j++ ) {
      tmp[j] = random.uniform( -1.0, 1.0 );
      noise[j] = random.normal( noiseStandardDeviation );
    }
            //project random point onto the sphere and add noise
    tmp.normalize();
    randomPoint = sphereCenter + noise + tmp*sphereRadius;
    data.push_back( randomPoint );
  }
           //generate outliers (via rejection)
  for( i=0; i<numOutliers; i++ ) {
    for( j=0; j<dimension; j++ ) {
      randomPoint[j] = random.uniform( -coordinateMax, coordinateMax );      
    }
    tmp = randomPoint - sphereCenter;
    if( tmp.l2Norm() >= outlierDistance )
      data.push_back( randomPoint );
    else
      i--;
  }
}

/*****************************************************************************/

template<unsigned int dimension>
void saveOIVFile(std::string &outputFileName,
                 std::vector< lsqrRecipes::Point<double,dimension> > &data, 
                 std::vector<double> &estimatedSphereParameters,
                 lsqrRecipes::SphereParametersEstimator<dimension> 
                   & parameterEstimator)
{
  if(dimension != 3)
    return;

  double dataSphereRadius = 50.0;
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

             //estimated sphere is gray
  std::string sphereMaterial = "\tMaterial {\n";
  sphereMaterial+="\t\tambientColor 0.5 0.5 0.5\n";
  sphereMaterial+="\t\tdiffuseColor 0.2 0.2 0.2\n";
  sphereMaterial+="\t\tspecularColor 0.8 0.8 0.8\n";
  sphereMaterial+="\ttransparency   0.4\n";
  sphereMaterial+="\t}\n";

  std::ofstream out(outputFileName.c_str());
  if(!out.is_open())
   return;

  out<<"#Inventor V2.1 ascii\n\n";
               //go over all the points and write 
               //the outliers and inliers to the scene graph
  for(unsigned int i=0; i<data.size(); i++) {
    out<<"Separator {\n";
    if(parameterEstimator.agree(estimatedSphereParameters, data[i])) {
      out<<inlierMaterial;     
    }
    else
      out<<outlierMaterial;
    out<<"\tTransform {\n";
    out<<"\t\ttranslation "<<(data[i])[0]<<" "<<(data[i])[1]<<" "<<(data[i])[2]<<"\n";
    out<<"\t}\n";
    out<<"\tSphere {\n";
    out<<"\t\tradius  "<<dataSphereRadius<<"\n";
    out<<"\t}\n";
    out<<"}\n";
  }
       //write estimated sphere to file  
  out<<"Separator {\n";
  out<<sphereMaterial;
  out<<"\tTransform {\n";
  out<<"\t\ttranslation "<<estimatedSphereParameters[0]<<" "<<estimatedSphereParameters[1]<<" "<<estimatedSphereParameters[2]<<"\n";
  out<<"\t}\n";
  out<<"\tSphere {\n";
  out<<"\t\tradius  "<<estimatedSphereParameters[3]<<"\n";
  out<<"\t}\n";
  out<<"}\n";
  
  out.close();
}
