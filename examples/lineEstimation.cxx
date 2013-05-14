#include <fstream>
#include "Point.h"
#include "Point3D.h"
#include "Vector.h"
#include "Vector3D.h"
#include "RANSAC.h"
#include "LineParametersEstimator.h"
#include "RandomNumberGenerator.h"

/**
 * Generate points on a kD line with additive zero mean Gaussian noise and
 * outliers.
 * @param numInliers How many points are inliers.
 * @param numOutliers How many points are outliers.
 * @param outlierDistance Threshold defining outliers, points that are farther
 *                        than this distance from the line.
 * @param data The points are added to the end of this vector.
 * @param lineParameters [n,a], line direction and point on line. Line is 
 *                        defined as the set of points p(t)=a+tn, where t is a
 *                        scalar in [-inf,inf]. 
 */
template<unsigned int dimension>
void generateData( unsigned int numInliers, unsigned int numOutliers, 
                   double outlierDistance,
                   std::vector< lsqrRecipes::Point<double,dimension> > &data,
                   std::vector<double> &lineParameters );

/**
 * This function only works in 3D.
 *
 * Save an open inventor ASCII scene file showing the points identified as
 * inliers/outliers and the estimated line (segment). The format is compatible 
 * with the viewers available in the open source coin3D toolkit (www.coin3d.org).
 *
 * @param outputFileName Write the scene to this file.
 * @param data The 3D points.
 * @param estimatedLineParameters [n,a], line direction and point on line. Line 
 *                                 is defined as the set of points p(t)=a+tn, 
 *                                 where t is a scalar in [-inf,inf].
 * @param parameterEstimator The line parameter estimator whose agree() method
 *                           is used to identify inliers.
 */
template<unsigned int dimension>
void saveOIVFile( std::string &outputFileName,
                  std::vector< lsqrRecipes::Point<double,dimension> > &data, 
                  std::vector<double> &estimatedLineParameters,
                  lsqrRecipes::LineParametersEstimator<dimension> 
                  &parameterEstimator );
/**
 * Given the hard coded dimension, and number of outliers and inliers generate
 * a random dataset accordingly. Then estimate the line parameter values
 * using a least squares estimate and the RANSAC algorithm. Compare the results
 * to the known line. Code is written for kD data except the 
 * visualization which is limited to 3D. If DIMENSION is set to three, two
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
  std::string leastSquaresOutputFileName = "leastSquaresLineEstimation.iv";
  std::string ransacOutputFileName = "RANSACLineEstimation.iv";

  std::vector< lsqrRecipes::Point<double,DIMENSION> > data;
  std::vector<double> trueLineParameters, lineParameters;
  lsqrRecipes::Point<double,DIMENSION> truePointOnLine, estimatedPointOnLine;
  lsqrRecipes::Vector<double,DIMENSION> trueLineDirection, tmp;

  double outlierDistance = 20.0;
  unsigned int i;  
  double dotProduct, distance;

  generateData<DIMENSION>(INLIERS, OUTLIERS, outlierDistance,
                          data, trueLineParameters);
  
  for( i=0; i<DIMENSION; i++ ) {
    trueLineDirection[i] = trueLineParameters[i];
    truePointOnLine[i] = trueLineParameters[DIMENSION+i];
  }
  std::cout<<"Known line parameters [n,a]\n\t [ ";
  for( i=0; i<(2*DIMENSION-1); i++ )
    std::cout<<trueLineParameters[i]<<", ";
  std::cout<<trueLineParameters[i]<<"]\n\n";


               //create and initialize the parameter estimator
                          //maximal distance from line is 2sigma of the Gaussian
                          //noise added to the point data in the generateData
                          //function below
  double maximalDistanceFromLine = 1;
  lsqrRecipes::LineParametersEstimator<DIMENSION> lineEstimator(maximalDistanceFromLine);
  lineEstimator.leastSquaresEstimate( data, lineParameters );
  if( lineParameters.empty() )
    std::cout<<"Least squares estimate failed, degenerate configuration?\n";
  else
  {
	  std::cout<<"Least squares line parameters: [n,a]\n\t [ ";
    for( i=0; i<(2*DIMENSION-1); i++ )
      std::cout<<lineParameters[i]<<", ";
    std::cout<<lineParameters[i]<<"]\n\n";

    for( i=0; i<DIMENSION; i++ )     
      estimatedPointOnLine[i] = lineParameters[DIMENSION+i];
      
                //cos(theta), theta is the angle between the two line directions
    dotProduct = 0.0;
    for( i=0; i<DIMENSION; i++ )
      dotProduct+= lineParameters[i]*trueLineParameters[i];
	  std::cout<<"\tDot product of real and computed line directions [+-1=correct]: ";
    std::cout<<dotProduct<<"\n";
              //distance between known line and estimated point on line  
    tmp = estimatedPointOnLine - truePointOnLine;
    distance = (tmp - tmp*trueLineDirection*trueLineDirection).l2Norm();    
	  std::cout<<"\tCheck if computed point is on known line [0=correct]: ";
    std::cout<<distance<<"\n\n";
             //save scene file (works only in 3D)
    saveOIVFile<DIMENSION>( leastSquaresOutputFileName, data, lineParameters, lineEstimator );
  }

                          //estimate using the RANSAC algorithm
  double desiredProbabilityForNoOutliers = 0.999;
  double percentageOfDataUsed;
  percentageOfDataUsed = 
    lsqrRecipes::RANSAC< lsqrRecipes::Point<double, DIMENSION>, double>::compute(lineParameters, 
                                                      &lineEstimator, 
                                                      data, 
                                                      desiredProbabilityForNoOutliers);
  if( lineParameters.empty() )
    std::cout<<"RANSAC estimate failed, degenerate configuration?\n";
  else
  {
	  std::cout<<"RANSAC line parameters: [n,a]\n\t [ ";
    for( i=0; i<(2*DIMENSION-1); i++ )
      std::cout<<lineParameters[i]<<", ";
    std::cout<<lineParameters[i]<<"]\n\n";

    for( i=0; i<DIMENSION; i++ )     
      estimatedPointOnLine[i] = lineParameters[DIMENSION+i];

                //cos(theta), theta is the angle between the two line directions
    dotProduct = 0.0;
    for( i=0; i<DIMENSION; i++ )
      dotProduct+= lineParameters[i]*trueLineParameters[i];
	  std::cout<<"\tDot product of real and computed line directions [+-1=correct]: ";
    std::cout<<dotProduct<<"\n";
              //distance between known line and estimated point on line  
    tmp = estimatedPointOnLine - truePointOnLine;
    distance = (tmp - tmp*trueLineDirection*trueLineDirection).l2Norm();    
	  std::cout<<"\tCheck if computed point is on known line [0=correct]: ";
    std::cout<<distance<<"\n\n";
    std::cout<<"\tPercentage of points which were used for final estimate: ";
    std::cout<<percentageOfDataUsed<<"\n\n";

             //save scene file (works only in 3D)
    saveOIVFile<DIMENSION>( ransacOutputFileName, data, lineParameters, lineEstimator );
  }
  return EXIT_SUCCESS;
}

/*****************************************************************************/

template<unsigned int dimension>
void generateData( unsigned int numInliers, unsigned int numOutliers, 
                   double outlierDistance, 
                   std::vector< lsqrRecipes::Point<double,dimension> > &data,
                   std::vector<double> &lineParameters )
{
  lsqrRecipes::Vector<double, dimension> lineDirection, tmp;
  lsqrRecipes::Point<double, dimension> pointOnLine, randomPoint;
  double noiseStandardDeviation = 0.5; //noise standard deviation
  double coordinateMax = 1000.0;
  double t;
  unsigned int i, j;

  lsqrRecipes::RandomNumberGenerator random;

  lineParameters.clear();
         //generate points on random line
  for( i=0; i<dimension; i++ ) {
    lineDirection[i] = random.uniform();
    pointOnLine[i] = random.uniform( -coordinateMax, coordinateMax );
  }
  lineDirection.normalize();
  for( i=0; i<dimension; i++ ) 
    lineParameters.push_back( lineDirection[i] ); 
  for( i=0; i<dimension; i++ ) 
    lineParameters.push_back( pointOnLine[i] );

               //generate inliers
  for( i=0; i<numInliers; i++ ) {
    t= random.uniform( -coordinateMax, coordinateMax );
    for( j=0; j<dimension; j++ )
      randomPoint[j] = pointOnLine[j] + t*lineDirection[j] + 
                       random.normal( noiseStandardDeviation );
    data.push_back( randomPoint );
  }
           //generate outliers (via rejection)
  for(i=0; i<numOutliers; i++) {
    for(j=0; j<dimension; j++) {
      randomPoint[j] = random.uniform( -coordinateMax, coordinateMax );      
    }
          //check that the distance of the random point is larger than the 
          //outlier distance
    tmp = randomPoint - pointOnLine;
    if( (tmp - tmp*lineDirection*lineDirection).l2Norm() >= outlierDistance )    
      data.push_back( randomPoint );
    else
      i--;
  }
}

/*****************************************************************************/

template<unsigned int dimension>
void saveOIVFile( std::string &outputFileName,
                  std::vector< lsqrRecipes::Point<double,dimension> > &data, 
                  std::vector<double> &estimatedLineParameters,
                  lsqrRecipes::LineParametersEstimator<dimension> 
                  &parameterEstimator )
{
  if( dimension != 3 )
    return;

  double sphereRadius = 10.0;
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

             //line has same color as inliers
  std::string lineMaterial = inlierMaterial;

  std::ofstream out( outputFileName.c_str() );
  if(!out.is_open())
   return;

  out<<"#Inventor V2.1 ascii\n\n";
               //go over all the points, find the coordinate ranges and write 
               //the outliers and inliers to the scene graph
  double minX, maxX, minY, maxY, minZ, maxZ;
  minX = maxX = data[0][0];
  minY = maxY = data[0][1];
  minZ = maxZ = data[0][2];
  unsigned int inlierIndex = 0;
  for( unsigned int i=0; i<data.size(); i++ ) {
    out<<"Separator {\n";
    if( parameterEstimator.agree( estimatedLineParameters, data[i] ) ) {
      out<<inlierMaterial;
      inlierIndex = i;
    }
    else
    out<<outlierMaterial;
    out<<"\tTransform {\n";
    out<<"\t\ttranslation "<<(data[i])[0]<<" "<<(data[i])[1]<<" "<<(data[i])[2]<<"\n";
    out<<"\t}\n";
    out<<"\tSphere {\n";
    out<<"\t\tradius  "<<sphereRadius<<"\n";
    out<<"\t}\n";
    out<<"}\n";
    if( data[i][0]<minX )
      minX = data[i][0];
    else if( data[i][0]>maxX )
      maxX = data[i][0];
    if( data[i][1]<minY )
      minY = data[i][1];
    else if( data[i][1]>maxY )
      maxY = data[i][1];
    if( data[i][2]<minZ )
      minZ = data[i][2];
    else if( data[i][2]>maxZ )
      maxZ = data[i][2];
  }
       //create line segment representing intersection of estimated line with
       //data's axis aligned bounding box, and write to file
      //the intersection point of of a line and a plane is given by:
      // N^T[a+tn-A]=0, where 'N' is the plane normal, 'A' a point on the plane,
      //'n' the line "direction", 'a' a point on the line, and 't' a scalar.
      //In our case we are working with an axis aligned box so all plane normals 
      //are of the form <1,0,0>. We can thus use the Liang-Barsky line clipping
      //algorithm in 3D
       
      //some sanity checks, just to make sure the line does indeed intersect
      //the bounding box
                       //line is parallel to the YZ plane
  if( fabs( estimatedLineParameters[0] ) < lsqrRecipes::EPS )  
		if( estimatedLineParameters[3]<minX || estimatedLineParameters[3] > maxX )
			throw std::exception();
                       //line is parallel to the XZ plane
  if( fabs( estimatedLineParameters[1] ) < lsqrRecipes::EPS )  
		if( estimatedLineParameters[4]<minY || estimatedLineParameters[4] > maxY )
			throw std::exception();
                       //line is parallel to the XY plane
  if( fabs( estimatedLineParameters[2] ) < lsqrRecipes::EPS )  
		if( estimatedLineParameters[5]<minZ || estimatedLineParameters[5] > maxZ )
			throw std::exception();
       
     //the line is "entering" the box if the dot product between the face normal
     //and the line direction is negative, if it is positive the line is exiting
     //the box. The parameter for the line intersection with the plane is given by:
     // t = (A_i - a_i)/n_i . 
  double tin, tout;
  int tEnteringNum, tLeavingNum,i; 
                 //arrays not initialized, which is fine, we keep tabs on the 
                 //number of valid entries
  double tEntering[3],tLeaving[3]; 
                                   
  tEnteringNum = tLeavingNum = 0;

                   //clip to right and left faces
  if( estimatedLineParameters[0] < 0 )
  {
		tEntering[tEnteringNum++] = (maxX-estimatedLineParameters[3])/estimatedLineParameters[0];
    tLeaving[tLeavingNum++] = (minX-estimatedLineParameters[3])/estimatedLineParameters[0];
  }
	else if( estimatedLineParameters[0] > 0 )
  {
    tEntering[tEnteringNum++] = (minX-estimatedLineParameters[3])/estimatedLineParameters[0];
		tLeaving[tLeavingNum++] = (maxX-estimatedLineParameters[3])/estimatedLineParameters[0];
  }

                   //clip to top and bottom faces
  if( estimatedLineParameters[1] < 0 )
  {
		tEntering[tEnteringNum++] = (maxY-estimatedLineParameters[4])/estimatedLineParameters[1];
    tLeaving[tLeavingNum++] = (minY-estimatedLineParameters[4])/estimatedLineParameters[1];
  }
	else if( estimatedLineParameters[1] > 0 )
  {
    tEntering[tEnteringNum++] = (minY-estimatedLineParameters[4])/estimatedLineParameters[1];
		tLeaving[tLeavingNum++] = (maxY-estimatedLineParameters[4])/estimatedLineParameters[1];
  }

                   //clip to front and back faces
  if( estimatedLineParameters[2] < 0 )
  {
		tEntering[tEnteringNum++] = (maxZ-estimatedLineParameters[5])/estimatedLineParameters[2];
    tLeaving[tLeavingNum++] = (minZ-estimatedLineParameters[5])/estimatedLineParameters[2];
  }
	else if( estimatedLineParameters[2] > 0 )
  {
    tEntering[tEnteringNum++] = (minZ-estimatedLineParameters[5])/estimatedLineParameters[2];
		tLeaving[tLeavingNum++] = (maxZ-estimatedLineParameters[5])/estimatedLineParameters[2];
  }
           //find the endpoints of the line segment intersecting the box
	tin = tEntering[0];
  for( i=0; i<tEnteringNum ; i++ ) {
	  if( tEntering[i] > tin ) {
		  tin = tEntering[i];
	  }
  }
  tout = tLeaving[0];
  for( i=0; i<tLeavingNum ; i++ ) {
	  if( tLeaving[i] < tout ) {
      tout = tLeaving[i];
	  }
  }
  
  out<<"Separator {\n";
  out<<lineMaterial;
  out<<"\tCoordinate3 {\n";
  out<<"\t\tpoint [";
  out<<estimatedLineParameters[3] + tin*estimatedLineParameters[0]<<" ";
  out<<estimatedLineParameters[4] + tin*estimatedLineParameters[1]<<" ";
  out<<estimatedLineParameters[5] + tin*estimatedLineParameters[2]<<", ";
  out<<estimatedLineParameters[3] + tout*estimatedLineParameters[0]<<" ";
  out<<estimatedLineParameters[4] + tout*estimatedLineParameters[1]<<" ";
  out<<estimatedLineParameters[5] + tout*estimatedLineParameters[2]<<"]\n";
  out<<"\t}\n";
  out<<"\tLineSet {\n";
  out<<"\t\tnumVertices [2]\n";
  out<<"\t}\n";
  out<<"}\n";

  out.close();
}
