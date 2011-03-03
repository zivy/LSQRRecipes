#include<vector>
#include <fstream>
#include <stdlib.h>
#include <vnl/vnl_cross.h>
#include "Point3D.h"
#include "RandomNumberGenerator.h"
#include "PlaneParametersEstimator.h"

using namespace lsqrRecipes;

/**
 * Generate points on a plane with additive Gaussian noise, except the first
 * three points.
 * @param numPoints How many points, in total there will be numPoints+3.
 * @param data The points are added to the end of this vector. The first three 
 *             points define the plane the rest numPoints points are 
 *             approximatly on the plane.
 */
void generatePlaneData(unsigned int numPoints, 
                       std::vector< Point3D > &data);
/**
 * Save the points to file.
 */
void saveData(const std::string &fileName, 
              std::vector< Point3D > &data);

/**
 * Load the points from file. It is assumed that the first three points
 * have no noise component, the rest do have additive Gaussian noise (see 
 * GeneratePlaneData()).
 */
void loadData(const std::string &fileName, 
              std::vector< Point3D > &data);

/*
 * Test the plane estimator's methods. We arbitrarily choose to work with 3D
 * planes even though the tested code works for any dimension.
 */
int main(int argc, char *argv[])
{
  const int NUM_SAMPLES = 20; //number of points sampled on plane
  std::vector< Point3D > pointData, minPointData;
  Point3D v1, v2, v3; 
  bool succeedExact, succeedLeastSquares;

  if(argc == 1) {
    generatePlaneData( NUM_SAMPLES, pointData );
    //saveData( std::string("planeTestData.txt"), pointData );
  }
  else
  {
    try {
      loadData( std::string( argv[1] ), pointData );
    } 
    catch( std::exception & ) {
      std::cerr<<"Failed to load input file.\n";
      return EXIT_FAILURE;
    }
  }
  minPointData.push_back(pointData[0]);
  minPointData.push_back(pointData[1]);
  minPointData.push_back(pointData[2]);
  
  

          //2. Test the code. Compare the known plane parameters to 
          //  (a) their exact estimate using three points; and 
          //  (b) their least squares estimate.

	std::vector<double> planeParameters;
  double maxDistanceToPlane = 0.5;
  PlaneParametersEstimator<3> ppEstimator(maxDistanceToPlane);

  double tmp;

	            //The known plane parameters
  v1 = pointData[0];
  v2 = pointData[1];
  v3 = pointData[2];
  vnl_vector<double> vec1(3), vec2(3), planeNormal(3);
  vec1[0] = v2[0] - v1[0];
  vec1[1] = v2[1] - v1[1];
  vec1[2] = v2[2] - v1[2];
  vec2[0] = v3[0] - v1[0];
  vec2[1] = v3[1] - v1[1];
  vec2[2] = v3[2] - v1[2];
  planeNormal = vnl_cross_3d(vec1, vec2);
  planeNormal.normalize();

	std::cout<<"Known plane parameters [nx,ny,nz,ax,ay,az]:\n\t [ ";
  std::cout<<planeNormal[0]<<", "<<planeNormal[1]<<", "<<planeNormal[2];
  std::cout<<", "<<v1[0]<<", "<<v1[1]<<", "<<v1[2]<<"]\n\n";

                  //check that the Agree method works
  planeParameters.push_back(planeNormal[0]);
  planeParameters.push_back(planeNormal[1]);
  planeParameters.push_back(planeNormal[2]);
  planeParameters.push_back(v1[0]);
  planeParameters.push_back(v1[1]);
  planeParameters.push_back(v1[2]);

  Point3D pointOnPlane = pointData[0];
  Point3D pointOffPlane;
  pointOffPlane[0] = v1[0] + 2 * maxDistanceToPlane * planeParameters[0];
  pointOffPlane[1] = v1[1] + 2 * maxDistanceToPlane * planeParameters[1];
  pointOffPlane[2] = v1[2] + 2 * maxDistanceToPlane * planeParameters[2];

  if( !ppEstimator.agree(planeParameters, pointOnPlane) ||
       ppEstimator.agree(planeParameters, pointOffPlane) )
    return EXIT_FAILURE;
  


             //compute an estimate using three points
  ppEstimator.estimate(minPointData, planeParameters);
	if(planeParameters.size() == 0) {
		std::cout<<"Plane going through three points [nx,ny,nz,ax,ay,az]: DEGENERATE CONFIGURATION\n\n";
    succeedExact = true;
  }
	else {
		std::cout<<"Plane going through three points [nx,ny,nz,ax,ay,az], no noise:\n\t [ ";
    std::cout<<planeParameters[0]<<", "<<planeParameters[1]<<", ";
		std::cout<<planeParameters[2]<<", "<<planeParameters[3]<<", ";
    std::cout<<planeParameters[4]<<", "<<planeParameters[5]<<" ]\n";		
    tmp = planeParameters[0]*planeNormal[0] + 
          planeParameters[1]*planeNormal[1] + 
          planeParameters[2]*planeNormal[2];
                   //angle between normals is less than 5 degrees
    succeedExact = fabs( tmp )  > 0.99619469809174553229501040247389;
		std::cout<<"\tDot product of known and computed plane normals[+-1=correct]: "<<tmp<<"\n";
    tmp = (planeParameters[3] - v1[0])*planeNormal[0] + 
          (planeParameters[4] - v1[1])*planeNormal[1] +
          (planeParameters[5] - v1[2])*planeNormal[2];
    succeedExact = succeedExact && tmp<maxDistanceToPlane;
		std::cout<<"\tTest if computed point is on known plane [0=correct]: "<<tmp<<"\n\n";
		planeParameters.clear();
	}
             //Least squares estimate of plane parameters
	ppEstimator.leastSquaresEstimate(pointData, planeParameters);
	if(planeParameters.size() == 0) {
		std::cout<<"Least squares plane parameters [nx,ny,nz,ax,ay,az]: DEGENERATE CONFIGURATION\n\n";
    succeedLeastSquares = true;
  }
	else {
		std::cout<<"Least squares plane parameters [nx,ny,nz,ax,ay,az]:\n\t [ ";
    std::cout<<planeParameters[0]<<", "<<planeParameters[1]<<", ";
		std::cout<<planeParameters[2]<<", "<<planeParameters[3]<<", ";
    std::cout<<planeParameters[4]<<", "<<planeParameters[5]<<" ]\n";		
    tmp = planeParameters[0]*planeNormal[0] + 
          planeParameters[1]*planeNormal[1] + 
          planeParameters[2]*planeNormal[2];
                   //angle between normals is less than 5 degrees
    succeedLeastSquares = fabs( tmp ) > 0.99619469809174553229501040247389;
		std::cout<<"\tDot product of known and computed plane normals[+-1=correct]: "<<tmp<<"\n";
    tmp = (planeParameters[3] - v1[0])*planeNormal[0] + 
          (planeParameters[4] - v1[1])*planeNormal[1] +
          (planeParameters[5] - v1[2])*planeNormal[2];
    succeedLeastSquares = succeedLeastSquares && tmp<maxDistanceToPlane;
		std::cout<<"\tTest if computed point is on known plane [0=correct]: "<<tmp<<"\n\n";
		planeParameters.clear();
	}
  if(succeedExact && succeedLeastSquares)
    return EXIT_SUCCESS;
  return EXIT_FAILURE;
}
/*****************************************************************************/
void generatePlaneData(unsigned int numPoints, 
                       std::vector< Point3D > &data)
{
	
  Point3D pnt;
  double v1[3], v2[3], v3[3];
  double bounds = 1000.0;

  double noiseSigma = 1.0; //noise is ~N(0,noiseSigma)

	      //1.Create data with noise: randomly select three points in
        //  [-bounds,bounds]x[-bounds,bounds]x[-bounds,bounds], then  
	      //  generate points on the plane defined by these points.
	      //  For each point sampled on the plane add random noise.

	RandomNumberGenerator random;

              //create random plane
  double threshold = 0.17453292519943295769236907684886; //sin(10^o)
  bool done = false;
  while(!done) {
               //three random points
	  v1[0] = random.uniform(-bounds, bounds);
    v1[1] = random.uniform(-bounds, bounds);
    v1[2] = random.uniform(-bounds, bounds);
	  v2[0] = random.uniform(-bounds, bounds);
    v2[1] = random.uniform(-bounds, bounds);
    v2[2] = random.uniform(-bounds, bounds);
	  v3[0] = random.uniform(-bounds, bounds);
    v3[1] = random.uniform(-bounds, bounds);
    v3[2] = random.uniform(-bounds, bounds);
                //compute plane normal and check for degenerate point 
                //configuration
    vnl_vector<double> vec1(3), vec2(3), planeNormal(3);
    vec1[0] = v2[0] - v1[0];
    vec1[1] = v2[1] - v1[1];
    vec1[2] = v2[2] - v1[2];
    vec2[0] = v3[0] - v1[0];
    vec2[1] = v3[1] - v1[1];
    vec2[2] = v3[2] - v1[2];
    planeNormal = vnl_cross_3d(vec1.normalize(), vec2.normalize());
    if(planeNormal.magnitude() > threshold)
      done = true;
  }
           //push the three points without adding noise
  data.push_back(Point3D(v1));
  data.push_back(Point3D(v2));
  data.push_back(Point3D(v3));

             //randomly generate points in the triangle defined by the above 
             //three points using barycentric coordinates, the number of 
             //generated points is numSamples            
  double w1,w2;             
  for(unsigned int i=0; i<numPoints; i++) {
    w1 = random.uniform(0.0, 1.0);
    w2 = random.uniform(0.0, 1.0);
           //we know that 0<=(w1+w2)<=2
    if(w1+w2>1) {
      w1 = 1-w1;
      w2 = 1-w2;
    }         
    pnt[0] = w1*v1[0] + w2*v2[0] + (1-w1-w2)*v3[0] + 
             random.normal(noiseSigma);
    pnt[1] = w1*v1[1] + w2*v2[1] + (1-w1-w2)*v3[1] + 
             random.normal(noiseSigma);
    pnt[2] = w1*v1[2] + w2*v2[2] + (1-w1-w2)*v3[2] + 
             random.normal(noiseSigma);
    data.push_back(pnt);
  } 
}
/*****************************************************************************/
void saveData(const std::string &fileName, 
              std::vector< Point3D > &data)
{
  std::ofstream out;
  out.open(fileName.c_str());
  if(out.fail()) 
    throw std::exception();
      
  for(unsigned int i=0; i<data.size(); i++) {
    out<<(data[i])[0]<<"\t"<<(data[i])[1]<<"\t"<<(data[i])[2]<<"\n";
  }
  out.close();
}
/*****************************************************************************/
void loadData(const std::string &fileName, 
              std::vector< Point3D > &data)
{
  Point3D pnt;

  std::ifstream in;
  in.open(fileName.c_str());
  if(in.fail()) 
    throw std::exception();
 
  while(!in.eof()) {
    in>>pnt[0]>>pnt[1]>>pnt[2];
    data.push_back(pnt);
  }
  in.close();
}
