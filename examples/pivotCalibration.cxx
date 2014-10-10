#include <fstream>
#include "RANSAC.h"
#include "PivotCalibrationParametersEstimator.h"


/**
 * Save an open inventor ASCII scene file showing the transformations identified 
 * as inliers/outliers and the estimated pivot calibration. The transformations 
 * are represented by coordinate systems where the "origin" is a green sphere 
 * for an inlier and a red one for an outlier. The pivot calibration is 
 * represented as a sphere with a radius of \|W^p\| and location \|{DRF}^p\|,
 * transparent red for the least squares estimate and green for the
 * RANSAC estimate.
 * The format is compatible with the viewers available in the open source 
 * coin3D toolkit (www.coin3d.org).
 *
 * @param outputFileName Write the scene to this file.
 * @param data The transformation frames.
 * @param leastSquaresEstimatedTranslations [{DRF}^t,W^t], results of pivot 
 *                                          calibration using the least squares method. 
 *                                          {DRF}^t is the translation from the DRF 
 *                                          to the pivoting point and W^t is
 *                                          the translation from the world/tracker coordinate
 *                                          system to the fixed pivoting point. 
 * @param ransacSquaresEstimatedTranslations [{DRF}^t,W^t], results of pivot 
 *                                          calibration using the RANSAC framework. 
 *                                          {DRF}^t is the translation from the DRF 
 *                                          to the pivoting point and W^t is
 *                                          the translation from the world/tracker coordinate
 *                                          system to the fixed pivoting point. 
 * @param pivotCalibration The pivot calibration estimator whose agree() method
 *                           is used to identify inliers.
 */
void saveOIVFile( std::string &outputFileName, 
                  std::vector< lsqrRecipes::Frame > &transformations, 
                  std::vector<double> &leastSquaresEstimatedTranslations,                   
                  std::vector<double> &ransacEstimatedTranslations,
                  lsqrRecipes::PivotCalibrationEstimator &pivotCalibration );

/**
 * Perform pivot calibration using transformations acquired with the Vicra optical
 * tracking system (Northern Digital Inc., Waterloo, Ontario, CA). The data consists
 * of 2/3 inliers and 1/3 outliers. We show that a least squares algebraic formulation
 * breaks down in the presence of outliers, but that by embedding it into the RANSAC
 * framework the correct results are obtained.
 * @author Ziv Yaniv
 */

int main( int argc, char *argv[] )
{
  if( argc != 2 ) {
    std::cerr<<"Unexpected number of command line parameters.\n";
    std::cerr<<"Usage: \n";
    std::cerr<<"\t"<<argv[0]<<" pivotCalibrationDataFileName\n";
    return EXIT_FAILURE;
  }

  double x, y, z, qx, qy, qz, qs;
  std::vector< lsqrRecipes::Frame > transformations;
  lsqrRecipes::Frame f;
  unsigned int i;
  const unsigned int PARAMETER_NUM = 6;
  std::string outputFileName = "PivotCalibration.iv";

               //load the data, each row in the file is expected to have the
               //following format x y z qx qy qz qs
  std::ifstream in;
  in.open( argv[1] );
  if( !in.is_open() )
    return EXIT_FAILURE;
  while( in>> x >> y >> z >> qx >> qy >> qz >> qs ) {
    f.setRotationQuaternion( qs, qx, qy, qz );
    f.setTranslation( x, y, z );   
    transformations.push_back( f );
  }
  in.close();

  if( transformations.empty() ) {
    std::cerr<<"Failed to load data from pivot calibration file ("<<argv[1]<<").\n";
    return EXIT_FAILURE;
  }

                        //create and initialize the parameter estimator
  double maxError = 1.0; //distance between the two points is at most 1mm
  std::vector<double> lsEstimatedTranslations, ransacEstimatedTranslations;
  lsqrRecipes::PivotCalibrationEstimator pivotCalibration( maxError );

                        //estimate using least squares
  pivotCalibration.leastSquaresEstimate( transformations, lsEstimatedTranslations );
  if( lsEstimatedTranslations.empty() ) {
    std::cout<<"Least squares estimate failed, degenerate configuration?\n";
    return EXIT_FAILURE;
  }
  else
  {
    std::cout<<"Least Squares estimate\n";
    std::cout<<"----------------------\n";
    std::cout<<"Estimated translations [DRF^t, W^t]:\n\t [ ";
    for( i=0; i<PARAMETER_NUM-1; i++ )
      std::cout<<lsEstimatedTranslations[i]<<", ";
    std::cout<<lsEstimatedTranslations[PARAMETER_NUM-1]<<"]\n\n";
  }
                          //estimate using the RANSAC algorithm
  double desiredProbabilityForNoOutliers = 0.999;
  double percentageOfDataUsed;
  percentageOfDataUsed = 
    lsqrRecipes::RANSAC< lsqrRecipes::Frame, double>::compute( ransacEstimatedTranslations, 
                                                               &pivotCalibration, 
                                                               transformations, 
                                                               desiredProbabilityForNoOutliers );
  if( ransacEstimatedTranslations.empty() ) {
    std::cout<<"RANSAC estimate failed, degenerate configuration?\n";
    return EXIT_FAILURE;
  }
  else
  {
    std::cout<<"RANSAC estimate\n";
    std::cout<<"---------------\n";
    std::cout<<"Estimated translations [DRF^t, W^t]:\n\t [ ";
    for( i=0; i<PARAMETER_NUM-1; i++ )
      std::cout<<ransacEstimatedTranslations[i]<<", ";
    std::cout<<ransacEstimatedTranslations[PARAMETER_NUM-1]<<"]\n\n";
  }

       //save scene file 
  saveOIVFile( outputFileName, transformations, lsEstimatedTranslations, 
               ransacEstimatedTranslations, pivotCalibration );

  return EXIT_SUCCESS;
}

/*****************************************************************************/

void saveOIVFile( std::string &outputFileName, 
                  std::vector< lsqrRecipes::Frame > &transformations, 
                  std::vector<double> &leastSquaresEstimatedTranslations,                   
                  std::vector<double> &ransacEstimatedTranslations,
                  lsqrRecipes::PivotCalibrationEstimator &pivotCalibration )
{
  const unsigned int SAMPLING_RATE = 10;

  std::ofstream out( outputFileName.c_str() );
  if( !out.is_open() )
   return;

  out<<"#Inventor V2.1 ascii\n\n";
  
           //the origin/sphere is green for an inlier
  std::string inlierString = "\tMaterial { diffuseColor 0 1 0 }\n";
  inlierString+="\tSphere {\n";
  inlierString+="\t\tradius 0.1\n";
  inlierString+="\t}\n";
  inlierString+="\tMaterial { diffuseColor 0 1 0 }\n";
  inlierString+="\tDEF arrow Separator\n";
  inlierString+="\t{\n";
  inlierString+="\t\tTranslation { translation 0 0.4 0 }\n";
  inlierString+="\t\tCylinder { radius 0.05\n";
  inlierString+="\t\t           height 0.8 }\n";
  inlierString+="\t\tTranslation { translation 0 0.4 0 }\n";
  inlierString+="\t\tCone { bottomRadius 0.1\n";
  inlierString+="\t\t       height 0.1 }\n";
  inlierString+="\t}\n";
  inlierString+="\tRotationXYZ { axis Z angle -1.5708 }\n";
  inlierString+="\tMaterial { diffuseColor 1 0 0 }\n";
  inlierString+="\tUSE arrow\n";
  inlierString+="\tRotationXYZ { axis X angle 1.5708 }\n";
  inlierString+="\tMaterial { diffuseColor 0 0 1 }\n";
  inlierString+="\tUSE arrow\n";
  inlierString+="}\n";

           //the origin/sphere is red for an outlier
  std::string outlierString = "\tMaterial { diffuseColor 1 0 0 }\n";
  outlierString+="\tSphere {\n";
  outlierString+="\t\tradius 0.1\n";
  outlierString+="\t}\n";
  outlierString+="\tMaterial { diffuseColor 0 1 0 }\n";
  outlierString+="\tDEF arrow Separator\n";
  outlierString+="\t{\n";
  outlierString+="\t\tTranslation { translation 0 0.4 0 }\n";
  outlierString+="\t\tCylinder { radius 0.05\n";
  outlierString+="\t\t           height 0.8 }\n";
  outlierString+="\t\tTranslation { translation 0 0.4 0 }\n";
  outlierString+="\t\tCone { bottomRadius 0.1\n";
  outlierString+="\t\t       height 0.1 }\n";
  outlierString+="\t}\n";
  outlierString+="\tRotationXYZ { axis Z angle -1.5708 }\n";
  outlierString+="\tMaterial { diffuseColor 1 0 0 }\n";
  outlierString+="\tUSE arrow\n";
  outlierString+="\tRotationXYZ { axis X angle 1.5708 }\n";
  outlierString+="\tMaterial { diffuseColor 0 0 1 }\n";
  outlierString+="\tUSE arrow\n";
  outlierString+="}\n";

  unsigned int i,n = transformations.size();
  double angleAxis[4], t[3]; 

  for( i=0; i<n; i+=SAMPLING_RATE ) {
    transformations[i].getRotationAxisAngle( angleAxis );
    transformations[i].getTranslation( t );
    out<<"Separator {\n";
    out<<"\tTransform { \n";
    out<<"\t\tscaleFactor "<<10<<"\t"<<10<<"\t"<<10<<"\n"; 
    out<<"\t\ttranslation "<<t[0]<<"\t"<<t[1]<<"\t"<<t[2]<<"\n"; 
    out<<"\t\trotation  "<<angleAxis[1]<<"\t"<<angleAxis[2]<<"\t"<<angleAxis[3]<<"\t"<<angleAxis[0]<<"\n";
    out<<"\t}\n";

    if( pivotCalibration.agree( ransacEstimatedTranslations, transformations[i] ) )
      out<<inlierString;
    else
      out<<outlierString;
  }    
              //display sphere with a radius of \|W^t\| and location \|{DRF}^t\|
              //transparent red for the least squares estimate and green for the
              //RANSAC estimate.
  double sphereRadius = 
    sqrt( ransacEstimatedTranslations[0]*ransacEstimatedTranslations[0]+
          ransacEstimatedTranslations[1]*ransacEstimatedTranslations[1]+
          ransacEstimatedTranslations[2]*ransacEstimatedTranslations[2] );
  out<<"Separator {\n";
  out<<"\tMaterial {\n";
  out<<"\t\tambientColor 0.0 1.0 0.0\n";
  out<<"\t}\n";
  out<<"\tTransform {\n";
  out<<"\t\ttranslation "<<ransacEstimatedTranslations[3]<<" "<<ransacEstimatedTranslations[4]<<" "<<ransacEstimatedTranslations[5]<<"\n";
  out<<"\t}\n";
  out<<"\tSphere {\n";
  out<<"\t\tradius  "<<sphereRadius<<"\n";
  out<<"\t}\n";
  out<<"}\n";

  sphereRadius = 
    sqrt( leastSquaresEstimatedTranslations[0]*leastSquaresEstimatedTranslations[0]+
          leastSquaresEstimatedTranslations[1]*leastSquaresEstimatedTranslations[1]+
          leastSquaresEstimatedTranslations[2]*leastSquaresEstimatedTranslations[2] );
  out<<"Separator {\n";
  out<<"\tMaterial {\n";
  out<<"\t\tambientColor 1.0 0.0 0.0\n";
  out<<"\t\ttransparency 0.5\n";
  out<<"\t}\n";
  out<<"\tTransform {\n";
  out<<"\t\ttranslation "<<leastSquaresEstimatedTranslations[3]<<" "<<leastSquaresEstimatedTranslations[4]<<" "<<leastSquaresEstimatedTranslations[5]<<"\n";
  out<<"\t}\n";
  out<<"\tSphere {\n";
  out<<"\t\tradius  "<<sphereRadius<<"\n";
  out<<"\t}\n";
  out<<"}\n";

  out.close();
}
