#include<vector>
#include <fstream>
#include <stdlib.h>
#include <vnl/vnl_cross.h>
#include "Point2D.h"
#include "RandomNumberGenerator.h"
#include "LineParametersEstimator.h"
#include "Line2DParametersEstimator.h"


using namespace lsqrRecipes;

/**
 * Generate points on a line with additive Gaussian noise, except the first
 * two points.
 * @param numPoints How many points, in total there will be numPoints+2.
 * @param data The points are added to the end of this vector. The first two 
 *             points define the line the rest numPoints points are 
 *             approximately on the line.
 */
void generateLineData(unsigned int numPoints, 
                      std::vector< Point2D > &data);
/**
 * Save the points to file.
 */
void saveData(const std::string &fileName, 
              std::vector< Point2D > &data);

/**
 * Load the points from file. It is assumed that the first two points
 * have no noise component, the rest do have additive Gaussian noise (see 
 * GenerateLineData()).
 */
void loadData(const std::string &fileName, 
              std::vector< Point2D > &data);

/*
 * Test the line estimator's methods. We choose to work with 2D
 * lines so that we test the code that is kD (LineParametersEstimator) and the 
 * specific 2D (Line2DParametersEstimator) in the same unit test. These two
 * estimators use slightly different parameterizations for the line. The kD 
 * estimator uses [n,a] where 'n' is a direction parallel to the line and 'a'
 * is a point on the line. The 2D specific version uses [n,a] where 'n' is the
 * normal to the line and 'a' is a point on the line. Only in 2D is the line 
 * normal unique (up to sign change), thus is parameterization is less appropriate
 * for kD.
 *
 * The testing code does appear repetitious (read copy paste), but this was the 
 * easy way to deal with the two types of parameterization.
 */
int main(int argc, char *argv[])
{
  const int NUM_SAMPLES = 20; //number of points sampled on the line
  std::vector< Point2D > pointData, minPointData;
  Point2D v1, v2; 
  bool succeedExact, succeedLeastSquares;

  if(argc == 1) {
    generateLineData( NUM_SAMPLES, pointData );
    //saveData( std::string("lineTestData.txt"), pointData );
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
  
  
          //2. Test the code. Compare the known line parameters to 
          //  (a) their exact estimate using two points; and 
          //  (b) their least squares estimate.

	std::vector<double> lineParameters;
  double maxDistanceToLine = 0.5;
  LineParametersEstimator<2> lpEstimator(maxDistanceToLine);
  Line2DParametersEstimator l2dpEstimator(maxDistanceToLine);
  double tmp;

	            //The known line parameters
  v1 = pointData[0];
  v2 = pointData[1];  
  vnl_vector<double> lineDirection(2);
  lineDirection[0] = v2[0] - v1[0];
  lineDirection[1] = v2[1] - v1[1];
  lineDirection.normalize();

	std::cout<<"Known line parameters [nx,ny,ax,ay], n is line direction:\n\t [ ";
  std::cout<<lineDirection[0]<<", "<<lineDirection[1];
  std::cout<<", "<<v1[0]<<", "<<v1[1]<<"]\n\n";

                  //check that the agree() method works
  lineParameters.push_back(lineDirection[0]);
  lineParameters.push_back(lineDirection[1]);
  lineParameters.push_back(v1[0]);
  lineParameters.push_back(v1[1]);

  Point2D pointOnLine = pointData[0];
  Point2D pointOffLine;
  pointOffLine[0] = v1[0] - 2 * maxDistanceToLine * lineParameters[1];
  pointOffLine[1] = v1[1] + 2 * maxDistanceToLine * lineParameters[0];
  
  if( !lpEstimator.agree(lineParameters, pointOnLine) ||
       lpEstimator.agree(lineParameters, pointOffLine) )
    return EXIT_FAILURE;
           //now check the 2D line estimator which uses a different 
           //parameterization [n,a] with 'n' perpendicular to the line direction 
  lineParameters[0] = - lineDirection[1];
  lineParameters[1] = lineDirection[0];
  if( !l2dpEstimator.agree(lineParameters, pointOnLine) ||
       l2dpEstimator.agree(lineParameters, pointOffLine) )
    return EXIT_FAILURE;


             //compute an estimate using two points
  lpEstimator.estimate(minPointData, lineParameters);
	if(lineParameters.size() == 0) {
		std::cout<<"Line going through two points [nx,ny,ax,ay]: DEGENERATE CONFIGURATION\n\n";
    succeedExact = true;
  }
	else {
		std::cout<<"Line going through two points [nx,ny,ax,ay], no noise:\n\t [ ";
    std::cout<<lineParameters[0]<<", "<<lineParameters[1]<<", ";
		std::cout<<lineParameters[2]<<", "<<lineParameters[3]<<" ]\n";		
    tmp = lineParameters[0]*lineDirection[0] + 
          lineParameters[1]*lineDirection[1];
                   //angle between line directions is less than 5 degrees
    succeedExact = fabs( tmp )  > 0.99619469809174553229501040247389;
		std::cout<<"\tDot product of known and computed line directions[+-1=correct]: "<<tmp<<"\n";
             //the normal to the line with direction [x,y] is [-y,x]
    tmp = -lineDirection[1]*(lineParameters[2] - v1[0]) + 
          lineDirection[0]*(lineParameters[3] - v1[1]);
    succeedExact = succeedExact && tmp<maxDistanceToLine;
		std::cout<<"\tTest if computed point is on known line [0=correct]: "<<tmp<<"\n\n";
		lineParameters.clear();
  }
             //repeat with the 2D specific version (different parameterization)
  l2dpEstimator.estimate(minPointData, lineParameters);
	if(lineParameters.size() == 0) {
		std::cout<<"Line going through two points [nx,ny,ax,ay]: DEGENERATE CONFIGURATION\n\n";
    succeedExact = succeedExact && true;
  }
	else {
		std::cout<<"Line going through two points [nx,ny,ax,ay], no noise:\n\t [ ";
    std::cout<<-lineParameters[1]<<", "<<lineParameters[0]<<", ";
		std::cout<<lineParameters[2]<<", "<<lineParameters[3]<<" ]\n";		
    tmp = -lineParameters[1]*lineDirection[0] + 
          lineParameters[0]*lineDirection[1];
                   //angle between line directions is less than 5 degrees
    succeedExact = succeedExact &&
                   fabs( tmp )  > 0.99619469809174553229501040247389;
		std::cout<<"\tDot product of known and computed line directions[+-1=correct]: "<<tmp<<"\n";
             //the normal to the line with direction [x,y] is [-y,x]
    tmp = -lineDirection[1]*(lineParameters[2] - v1[0]) + 
          lineDirection[0]*(lineParameters[3] - v1[1]);
    succeedExact = succeedExact && tmp<maxDistanceToLine;
		std::cout<<"\tTest if computed point is on known line [0=correct]: "<<tmp<<"\n\n";
		lineParameters.clear();
  }

             //Least squares estimate of line parameters
	lpEstimator.leastSquaresEstimate(pointData, lineParameters);
	if(lineParameters.size() == 0) {
		std::cout<<"Least squares line parameters [nx,ny,ax,ay]: DEGENERATE CONFIGURATION\n\n";
    succeedLeastSquares = true;
  }
	else {
		std::cout<<"Least squares line parameters [nx,ny,ax,ay]:\n\t [ ";
    std::cout<<lineParameters[0]<<", "<<lineParameters[1]<<", ";
		std::cout<<lineParameters[2]<<", "<<lineParameters[3]<<" ]\n";		
    tmp = lineParameters[0]*lineDirection[0] + 
          lineParameters[1]*lineDirection[1];
                   //angle between line directions is less than 5 degrees
    succeedLeastSquares = fabs( tmp ) > 0.99619469809174553229501040247389;
		std::cout<<"\tDot product of known and computed line directions [+-1=correct]: "<<tmp<<"\n";
             //the normal to the line with direction [x,y] is [-y,x]
    tmp = -lineDirection[1]*(lineParameters[2] - v1[0]) + 
          lineDirection[0]*(lineParameters[3] - v1[1]);
    succeedLeastSquares = succeedLeastSquares && tmp<maxDistanceToLine;
		std::cout<<"\tTest if computed point is on known line [0=correct]: "<<tmp<<"\n\n";
		lineParameters.clear();
	}

             //repeat with the 2D specific version (different parameterization)
  l2dpEstimator.leastSquaresEstimate(pointData, lineParameters);
	if(lineParameters.size() == 0) {
		std::cout<<"Least squares line parameters [nx,ny,ax,ay]: DEGENERATE CONFIGURATION\n\n";
    succeedLeastSquares = succeedLeastSquares && true;
  }
	else {
		std::cout<<"Least squares line parameters [nx,ny,ax,ay]:\n\t [ ";
    std::cout<<-lineParameters[1]<<", "<<lineParameters[0]<<", ";
		std::cout<<lineParameters[2]<<", "<<lineParameters[3]<<" ]\n";		
    tmp = -lineParameters[1]*lineDirection[0] + 
          lineParameters[0]*lineDirection[1];
                   //angle between line directions is less than 5 degrees
    succeedLeastSquares = succeedLeastSquares &&
                          fabs( tmp ) > 0.99619469809174553229501040247389;
		std::cout<<"\tDot product of known and computed line directions [+-1=correct]: "<<tmp<<"\n";
             //the normal to the line with direction [x,y] is [-y,x]
    tmp = -lineDirection[1]*(lineParameters[2] - v1[0]) + 
          lineDirection[0]*(lineParameters[3] - v1[1]);
    succeedLeastSquares = succeedLeastSquares && tmp<maxDistanceToLine;
		std::cout<<"\tTest if computed point is on known line [0=correct]: "<<tmp<<"\n\n";
		lineParameters.clear();
	}

  if(succeedExact && succeedLeastSquares)
    return EXIT_SUCCESS;
  return EXIT_FAILURE;
}
/*****************************************************************************/
void generateLineData(unsigned int numPoints, 
                      std::vector< Point2D > &data)
{
	
  Point2D pnt;
  double a[2], n[2], p[2], t;
  double directionNorm, bounds = 1000.0;

  double noiseSigma = 1.0; //noise is ~N(0,noiseSigma)

	      //1.Create data with noise: randomly select a point 'a' in
        //  [-bounds,bounds]x[-bounds,bounds] and a direction 'n', then  
	      //  generate points on a line segment defined by a+tn where 't' is a 
        //  scalar in [-bounds, bounds].
	      //  For each point sampled on the line add random noise.

	RandomNumberGenerator random;

              //create random line               
	 a[0] = random.uniform(-bounds, bounds);
   a[1] = random.uniform(-bounds, bounds);
	 n[0] = random.uniform(-1.0, 1.0);
   n[1] = random.uniform(-1.0, 1.0);
   directionNorm = sqrt(n[0]*n[0] + n[1]*n[1]);
   n[0]/=directionNorm;
   n[1]/=directionNorm;

           //push two points without adding noise
  data.push_back(Point2D(a));
  t = random.uniform(-bounds, bounds);
  p[0] = a[0] + t*n[0];
  p[1] = a[1] + t*n[1];
  data.push_back(Point2D(p));

             //randomly generate points on the line defined above 
             //adding normally distributed noise to the coordinates, the number 
             //of generated points is numSamples            
  for(unsigned int i=0; i<numPoints; i++) {
    t = random.uniform(-bounds, bounds);
    p[0] = a[0] + t*n[0] + random.normal(noiseSigma);
    p[1] = a[1] + t*n[1] + random.normal(noiseSigma);
    data.push_back(Point2D(p));
  }
}
/*****************************************************************************/
void saveData(const std::string &fileName, 
              std::vector< Point2D > &data)
{
  std::ofstream out;
  out.open(fileName.c_str());
  if(out.fail()) 
    throw std::exception();
      
  for(unsigned int i=0; i<data.size(); i++) {
    out<<(data[i])[0]<<"\t"<<(data[i])[1]<<"\n";
  }
  out.close();
}
/*****************************************************************************/
void loadData(const std::string &fileName, 
              std::vector< Point2D > &data)
{
  Point2D pnt;

  std::ifstream in;
  in.open(fileName.c_str());
  if(in.fail()) 
    throw std::exception();
 
  while(!in.eof()) {
    in>>pnt[0]>>pnt[1];
    data.push_back(pnt);
  }
  in.close();
}
