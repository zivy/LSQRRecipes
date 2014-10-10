#include <fstream>
#include <stdlib.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_rank.h>
#include "RandomNumberGenerator.h"
#include "DenseLinearEquationSystemParametersEstimator.h"

using namespace lsqrRecipes;

/**
 * Use simulated data, clean data to test the exact estimate and data with
 * noise to test the least squares estimate.
 */
int simulatedDataTest();

/**
 * Use experimental data acquired during pivot calibration to test the least
 * squares method. This deals with actual noise and not noise artificially 
 * generated using a noise model we selected.
 */
int pivotCalibrationMatrixTest(const std::string &augmentedMatrixFileName);

int main(int argc, char *argv[])
{
  if(argc == 1)
    return simulatedDataTest();
  else if(argc == 2)
    return pivotCalibrationMatrixTest(argv[1]);
  else {
    std::cerr<<"Unexpected number of command line parameters.\n";
    std::cerr<<"Usage: \n";
    std::cerr<<"\t"<<argv[0]<<"\n";
    std::cerr<<"\t\t\tor\n";
    std::cerr<<"\t"<<argv[0]<<" augmentedMatrixFileName\n";
    return EXIT_FAILURE;
  }
}

/*****************************************************************************/

int simulatedDataTest()
{
  bool succeedExact, succeedLeastSquares;
  const unsigned int rowNum = 200;
  const unsigned int colNum = 5;
  unsigned int i, j;
  bool ok;

	vnl_matrix<double> A( rowNum, colNum ), ASquare( colNum, colNum );
  vnl_vector<double> b( rowNum ), bNoise( rowNum ), bSquare( colNum), x( colNum );

  RandomNumberGenerator random;

        //random invertible matrix
  ok = false;
  while (!ok) {
    for( i=0; i<colNum; i++ ) {
      for( j=0; j<colNum; j++ ) {
        ASquare(j,i) = random.uniform( -1.0, 1.0 );
      }
      x(i) = random.uniform( -1.0, 1.0 );
    }
    ok = vnl_rank( ASquare )== colNum;
  }
  bSquare = ASquare*x;

  std::vector< AugmentedRow<double, colNum> > augmentedRows;
  DenseLinearEquationSystemParametersEstimator<double,colNum>::getAugmentedRows( ASquare, bSquare, 
                               augmentedRows );
  
  double maxEquationError = 1.0e-10;
  std::vector<double> estimatedSolution;
  DenseLinearEquationSystemParametersEstimator<double, colNum> equationSolver( maxEquationError );

  equationSolver.estimate( augmentedRows, estimatedSolution ); 
  if(estimatedSolution.empty())
    return EXIT_FAILURE;

  std::cout<<"Exact solution with invertible matrix\n";
  std::cout<<"-------------------------------------\n";
  std::cout<<"Known solution [x_0,..., x_{n-1}]:\n\t [ ";
  for( i=0; i<colNum-1; i++ )
    std::cout<<x[i]<<", ";
  std::cout<<x[colNum-1]<<"]\n";

  std::cout<<"Estimated solution [x_0,..., x_{n-1}]:\n\t [ ";
  for( i=0; i<colNum-1; i++ )
    std::cout<<estimatedSolution[i]<<", ";
  std::cout<<estimatedSolution[colNum-1]<<"]\n\n";

                 //solution without any noise should have no error, so we use the
                 //tight bound defined by maxEquationError
  succeedExact = true;
  for( i=0; i<colNum; i++ )
    succeedExact = succeedExact && (fabs(estimatedSolution[i] - x[i])<maxEquationError);
  
  //random over-determined system, with rank colNum
  ok = false;
  while (!ok) {
    for( i=0; i<colNum; i++ ) {
      for( j=0; j<rowNum; j++ ) {
        A(j,i) = random.uniform( -1.0, 1.0 );
      }
      x(i) = random.uniform( -1.0, 1.0 );
    }
    ok = vnl_rank( A )== colNum;
  }
  b = A*x;

            //modify b by adding noise, where the noise is a specific percentile 
            //of the signal
  double noiseBounds, maxNoiseScale=0.05;

  for( i=0; i<rowNum; i++ ) {
    noiseBounds = fabs( maxNoiseScale*b(i) );
    bNoise(i) = random.uniform( -noiseBounds, noiseBounds );
  }
  b += bNoise;

  DenseLinearEquationSystemParametersEstimator<double,colNum>::getAugmentedRows(A, b, 
                               augmentedRows);
  maxEquationError = 0.1;
  equationSolver.setDelta(maxEquationError);
  equationSolver.leastSquaresEstimate(augmentedRows, estimatedSolution); 
  if(estimatedSolution.empty())
    return EXIT_FAILURE;

  std::cout<<"Least squares solution with overdetermined  matrix and noise\n";
  std::cout<<"------------------------------------------------------------\n";
  std::cout<<"Known solution [x_0,..., x_{n-1}]:\n\t [ ";
  for( i=0; i<colNum-1; i++ )
    std::cout<<x[i]<<", ";
  std::cout<<x[colNum-1]<<"]\n";

  std::cout<<"Estimated solution [x_0,..., x_{n-1}]:\n\t [ ";
  for( i=0; i<colNum-1; i++ )
    std::cout<<estimatedSolution[i]<<", ";
  std::cout<<estimatedSolution[colNum-1]<<"]\n\n";

                 //least squares solution with noise is expected to have errors in 
                 //so we use the loose bound re-defined by maxEquationError
  succeedLeastSquares = true;
  for( i=0; i<colNum; i++ ) {
    succeedLeastSquares = succeedLeastSquares && (fabs(estimatedSolution[i] - x[i])<maxEquationError);    
  }

  if(succeedExact && succeedLeastSquares)
    return EXIT_SUCCESS;
  return EXIT_FAILURE;
}

/*****************************************************************************/

int pivotCalibrationMatrixTest(const std::string &augmentedMatrixFileName)
{
  const unsigned int COL_NUM = 6;
  double data[COL_NUM+1];
  AugmentedRow<double, COL_NUM> row;
  std::vector< AugmentedRow<double, COL_NUM> > augmentedRows;
     //known solution to equation system
  double x[6] = { -1.777985584409468e+001, 1.111302171667757e+000, 
                  -1.568653413096010e+002, 1.469013927556186e+002, 
                  -6.296891425314718e+001,-1.042139650090033e+003};

    //load the data
  std::ifstream in;
  in.open(augmentedMatrixFileName.c_str());
  if(!in.is_open())
    return false;
  while(in>>data[0]>>data[1]>>data[2]>>data[3]>>data[4]>>data[5]>>data[6]) {
    row.set(data);    
    augmentedRows.push_back(row);
  }
  in.close();

  if(augmentedRows.empty()) {
    std::cerr<<"Failed to load augmented matrix file.\n";
    return EXIT_FAILURE;
  }

      //set the equation solver
  double maxEquationError = 0.5;
  std::vector<double> estimatedSolution;
  DenseLinearEquationSystemParametersEstimator<double, COL_NUM> equationSolver( maxEquationError );

  equationSolver.leastSquaresEstimate(augmentedRows, estimatedSolution);
  if(estimatedSolution.empty())
    return EXIT_FAILURE;

  unsigned int i;
  
  std::cout<<"Least squares solution experimental overdetermined  matrix\n";
  std::cout<<"----------------------------------------------------------\n";
  std::cout<<"Known solution [x_0,..., x_{n-1}]:\n\t [ ";
  for( i=0; i<COL_NUM-1; i++ )
    std::cout<<x[i]<<", ";
  std::cout<<x[COL_NUM-1]<<"]\n";

  std::cout<<"Estimated solution [x_0,..., x_{n-1}]:\n\t [ ";
  for( i=0; i<COL_NUM-1; i++ )
    std::cout<<estimatedSolution[i]<<", ";
  std::cout<<estimatedSolution[COL_NUM-1]<<"]\n\n";

  bool succeedLeastSquares = true;
  for( i=0; i<COL_NUM; i++ ) {
    succeedLeastSquares = succeedLeastSquares && (fabs(estimatedSolution[i] - x[i])<maxEquationError);    
  }
  return succeedLeastSquares ? EXIT_SUCCESS : EXIT_FAILURE;
}
