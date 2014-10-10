#ifndef _DENSE_LINEAR_EQUATION_SYSTEM_PARAMETERS_ESTIMATOR_TXX_ 
#define _DENSE_LINEAR_EQUATION_SYSTEM_PARAMETERS_ESTIMATOR_TXX_

#include<vnl/algo/vnl_matrix_inverse.h>
#include "Epsilon.h"

namespace lsqrRecipes {

template<class T, unsigned int n> 	
DenseLinearEquationSystemParametersEstimator<T,n>::DenseLinearEquationSystemParametersEstimator( T delta ) : 
  ParametersEstimator< AugmentedRow<T, n>, T>(n) 
{
  this->delta = delta;
}
/*****************************************************************************/
template<class T, unsigned int n> 	
void DenseLinearEquationSystemParametersEstimator<T,n>::estimate(
  std::vector< AugmentedRow<T, n> *> &data, 
  std::vector<T> &parameters )
{
  vnl_matrix<T> A(n,n);
  vnl_vector<T> x(n), b(n);
  T rowData[n], bValue;
  unsigned int i, numRows = static_cast<unsigned int>(data.size());  

  parameters.clear();
  	              //not enough data elements for computation
	if(numRows<this->minForEstimate)
		return;

              //set up the equation system
  for(i=0; i<numRows; i++) {
    data[i]->get(rowData, bValue);
    A.set_row(i, rowData);
    b[i] = bValue;
  }

  vnl_matrix_inverse<double> Ainv(A);
                     //explicitly zero out small singular values
                     //this is ugly as it exposes that the inverse is computed via SVD
  Ainv.zero_out_absolute(EPS);

  if(Ainv.rank()<n) //matrix is not invertible
    return;
  x = Ainv * b;
  for(i=0; i<n; i++) {
    parameters.push_back(x[i]);
  }  
}
/*****************************************************************************/
template<class T, unsigned int n> 	
void DenseLinearEquationSystemParametersEstimator<T,n>::estimate(
  std::vector< AugmentedRow<T, n> > &data, 
  std::vector<T> &parameters )
{
	std::vector< AugmentedRow<T, n> *> usedData;
	unsigned int dataSize = static_cast<unsigned int>(data.size());  
	for(unsigned int i=0; i<dataSize; i++)
		usedData.push_back(&(data[i]));
	estimate(usedData,parameters);
}
/*****************************************************************************/
template<class T, unsigned int n> 	
void DenseLinearEquationSystemParametersEstimator<T,n>::leastSquaresEstimate(
  std::vector< AugmentedRow<T, n> *> &data, 
  std::vector<T> &parameters )
{
  unsigned int i, numRows = static_cast<unsigned int>(data.size());  
  vnl_matrix<T> A(numRows,n);
  vnl_vector<T> x(n), b(numRows);
  T rowData[n], bValue;
  
  parameters.clear();
  	              //not enough data elements for computation
	if(numRows<this->minForEstimate)
		return;

              //set up the equation system
  for(i=0; i<numRows; i++) {
    data[i]->get(rowData, bValue);
    A.set_row(i, rowData);
    b[i] = bValue;
  }

  vnl_matrix_inverse<double> Ainv(A);
                     //explicitly zero out small singular values
                     //this is ugly as it exposes that the inverse is computed via SVD
  Ainv.zero_out_absolute(EPS);

  if(Ainv.rank()<n) //we do not have a unique solution 
    return;
  x = Ainv * b;
  for(i=0; i<n; i++) {
    parameters.push_back(x[i]);
  }  
}
/*****************************************************************************/
template<class T, unsigned int n> 	
void DenseLinearEquationSystemParametersEstimator<T,n>::leastSquaresEstimate(
  std::vector< AugmentedRow<T, n> > &data, 
  std::vector<T> &parameters )
{
	std::vector< AugmentedRow<T, n> *> usedData;
	unsigned int dataSize = static_cast<unsigned int>(data.size());  
	for(unsigned int i=0; i<dataSize; i++)
		usedData.push_back(&(data[i]));
	leastSquaresEstimate(usedData,parameters);
}
/*****************************************************************************/
template<class T, unsigned int n>
bool DenseLinearEquationSystemParametersEstimator<T,n>::agree(
  std::vector<T> &parameters, AugmentedRow<T, n> &data )
{
  T sum = static_cast<T>(0.0);
  for( unsigned int i=0; i<n; i++ )
    sum += data[i]*parameters[i];
  sum-=data[n];
  return fabs( sum ) < this->delta;
}
/*****************************************************************************/
template<class T, unsigned int n>
void DenseLinearEquationSystemParametersEstimator<T,n>::getAugmentedRows(
  vnl_matrix<T> &A, vnl_vector<T> &b, std::vector< AugmentedRow<T, n> > &rows )
{
  unsigned int rowNum = A.rows();
         //ensure the vector and matrix dimensions are consistent and that they
         //are consistent with the output
  if( b.size() != rowNum || A.cols() != n )
    throw std::exception();
  
  rows.resize( rowNum );  
  for( unsigned int i=0; i<rowNum; i++ ) {
    rows[i].set(A.get_row(i).data_block(), b(i));
  }

}


} //namespace lsqrRecipes

#endif //_DENSE_LINEAR_EQUATION_SYSTEM_PARAMETERS_ESTIMATOR_TXX_
