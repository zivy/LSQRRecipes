#include <fstream>
#include <vnl/vnl_rank.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include "RANSAC.h"
#include "RandomNumberGenerator.h"
#include "DenseLinearEquationSystemParametersEstimator.h"
#include "Vector.h"

/**
 * Create a random equation system, scale a number of the equations so that they
 * become outliers and estimate the solution. This system only includes outliers
 * and no noise.
 */
void simulatedExample();

/**
 * Estimate the solution of an equation system that includes both noise and 
 * outliers. This equation system solves the pivot calibration task, with data
 * acquired using an optically tracked pointer tool. We intentionally introduced
 * 30% outliers by arbitrarily rotation and translating the tool.
 */
void experimentalExample(char *inputFileName);


/**
 * Show how to use the RANSAC framework as a robust approach to solving a set of
 * linear equations.
 *
 * @author Ziv Yaniv 
 */
int main(int argc, char *argv[])
{
  simulatedExample();
  if(argc==2)
    experimentalExample(argv[1]);
}

/*****************************************************************************/

void simulatedExample()
{
  const unsigned int ROW_NUM = 200;
  const unsigned int COL_NUM = 5;
  double scaleFactor = 20.0;
  bool ok;
  unsigned int i,j;

  const unsigned int OUTLIER_NUM = 5;
                      // indexes to equations that will be outliers (for the fun
                      //of it use prime numbers)
  unsigned int outlierIndexes[OUTLIER_NUM] = {13, 37, 67, 137, 179}; 

  vnl_matrix<double> A(ROW_NUM, COL_NUM);
  vnl_vector<double> b(ROW_NUM), x(COL_NUM);
  lsqrRecipes::Vector<double, COL_NUM> tmp; 
  lsqrRecipes::RandomNumberGenerator random;

           //display two digits after the decimal point, settings are returned
           //to their original values at the end of the method
  std::ios::fmtflags previousFlags = 
    std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::streamsize previousPrecision = std::cout.precision(2);

  ok = false;
  while (!ok) {
    for( i=0; i<COL_NUM; i++ ) {
      for( j=0; j<ROW_NUM; j++ ) {
        A(j,i) = random.uniform( -100, 100 );
      }
      x(i) = random.uniform( -100, 100 );
    }
    ok = vnl_rank( A )== COL_NUM;
  }
  b = A*x;

       //create several outlying equations
  for(i=0; i<OUTLIER_NUM; i++)
    for(j=0; j<COL_NUM; j++ )
      A(outlierIndexes[i], j) *= scaleFactor;

       //translate the data into the format used by our estimator   
  std::vector< lsqrRecipes::AugmentedRow<double, COL_NUM> > augmentedRows;
  lsqrRecipes::DenseLinearEquationSystemParametersEstimator<double,COL_NUM>::getAugmentedRows( A, b, 
                               augmentedRows );

  std::cout<<"Simulated example (without noise, with outliers)\n";
  std::cout<<"------------------------------------------------\n";

  std::cout<<"Known solution [x_0,..., x_{n-1}]:\n\t [ ";
  for( i=0; i<COL_NUM-1; i++ )
    std::cout<<x[i]<<", ";
  std::cout<<x[COL_NUM-1]<<"]\n";

                        //create and initialize the parameter estimator
  double maximalEquationError = 0.2;
  std::vector<double> estimatedSolution;
  lsqrRecipes::DenseLinearEquationSystemParametersEstimator<double, COL_NUM> equationSystemSolver(maximalEquationError);

                     //estimate using a least squares approach
  equationSystemSolver.leastSquaresEstimate(augmentedRows, estimatedSolution);
  if( estimatedSolution.empty() )
    std::cout<<"Least squares estimate failed, degenerate configuration?\n";
  else {
    std::cout<<"Least squares estimated solution [x_0,..., x_{n-1}]:\n\t [ ";
    for(i=0; i<COL_NUM-1; i++) {
      std::cout<<estimatedSolution[i]<<", ";
      tmp[i] = estimatedSolution[i] - x[i];
    }
    std::cout<<estimatedSolution[COL_NUM-1]<<"]\n";
    tmp[i] = estimatedSolution[COL_NUM-1] - x[COL_NUM-1];
    std::cout<<"\t Distance between estimated and known solutions [0=correct]: ";
    std::cout<<tmp.l2Norm()<<"\n";
  } 

                          //estimate using the RANSAC algorithm
  double desiredProbabilityForNoOutliers = 0.999;
  double percentageOfDataUsed;
  percentageOfDataUsed = 
    lsqrRecipes::RANSAC< lsqrRecipes::AugmentedRow<double, COL_NUM>, double>::compute(estimatedSolution, 
                                                      &equationSystemSolver, 
                                                      augmentedRows, 
                                                      desiredProbabilityForNoOutliers);
  if( estimatedSolution.empty() )
    std::cout<<"RANSAC estimate failed, degenerate configuration?\n";
  else
  {
	  std::cout<<"RANSAC estimated solution [x_0,..., x_{n-1}]:\n\t [ ";
    for(i=0; i<COL_NUM-1; i++) {
      std::cout<<estimatedSolution[i]<<", ";
      tmp[i] = estimatedSolution[i] - x[i];
    }
    std::cout<<estimatedSolution[COL_NUM-1]<<"]\n";
    tmp[i] = estimatedSolution[COL_NUM-1] - x[COL_NUM-1];
    std::cout<<"\t Distance between estimated and known solutions [0=correct]: ";
    std::cout<<tmp.l2Norm()<<"\n\n";
  }
       //return the stream settings to their original values    
  std::cout.setf(previousFlags);
  std::cout.precision(previousPrecision);
}

/*****************************************************************************/

void experimentalExample(char *inputFileName)
{
  const unsigned int COL_NUM = 6;
  double data[COL_NUM+1];
  lsqrRecipes::AugmentedRow<double, COL_NUM> row;
  std::vector< lsqrRecipes::AugmentedRow<double, COL_NUM> > augmentedRows;
  unsigned int i;

  std::ifstream in;
  in.open(inputFileName);
  if(!in.is_open()) {
    std::cerr<<"Failed to open input file ("<<inputFileName<<").\n";
    return;
  }
  while(in>>data[0]>>data[1]>>data[2]>>data[3]>>data[4]>>data[5]>>data[6]) {
    row.set(data);    
    augmentedRows.push_back(row);
  }
  in.close();

  if(augmentedRows.empty()) {
    std::cerr<<"Failed to load augmented matrix file.\n";
    return;
  }

           //display two digits after the decimal point, settings are returned
           //to their original values at the end of the method
  std::ios::fmtflags previousFlags = 
    std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::streamsize previousPrecision = std::cout.precision(2);

  std::cout<<"Experimental example acquired for pivot calibration (noise + outliers)\n";
  std::cout<<"----------------------------------------------------------------------\n";

  std::cout<<"Correct and approximately known solution [x_0,..., x_{n-1}]:\n\t";
  std::cout<<"[ -17.0, 1.0, -157.0, 147.0, -63.0, -1042.0]\n";

                        //create and initialize the parameter estimator
  double maximalEquationError = sqrt(static_cast<double>(1.0/3.0));//0.3;
  std::vector<double> estimatedSolution;
  lsqrRecipes::DenseLinearEquationSystemParametersEstimator<double, COL_NUM> equationSystemSolver(maximalEquationError);

                     //estimate using a least squares approach
  equationSystemSolver.leastSquaresEstimate(augmentedRows, estimatedSolution);
  if( estimatedSolution.empty() )
    std::cout<<"Least squares estimate failed, degenerate configuration?\n";
  else {
    std::cout<<"Least squares estimated solution [x_0,..., x_{n-1}]:\n\t [ ";
    for(i=0; i<COL_NUM-1; i++) 
      std::cout<<estimatedSolution[i]<<", ";
    std::cout<<estimatedSolution[COL_NUM-1]<<"]\n";
  } 

                          //estimate using the RANSAC algorithm
  double desiredProbabilityForNoOutliers = 0.999;
  double percentageOfDataUsed;
  percentageOfDataUsed = 
    lsqrRecipes::RANSAC< lsqrRecipes::AugmentedRow<double, COL_NUM>, double>::compute(estimatedSolution, 
                                                      &equationSystemSolver, 
                                                      augmentedRows, 
                                                      desiredProbabilityForNoOutliers);
  if( estimatedSolution.empty() )
    std::cout<<"RANSAC estimate failed, degenerate configuration?\n";
  else
  {
	  std::cout<<"RANSAC estimated solution [x_0,..., x_{n-1}]:\n\t [ ";
    for(i=0; i<COL_NUM-1; i++) 
      std::cout<<estimatedSolution[i]<<", ";
    std::cout<<estimatedSolution[COL_NUM-1]<<"]\n";
  }

       //return the stream settings to their original values    
  std::cout.setf(previousFlags);
  std::cout.precision(previousPrecision);
}