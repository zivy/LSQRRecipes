#ifndef _LINE_PARAMETERS_ESTIMATOR_TXX_
#define _LINE_PARAMETERS_ESTIMATOR_TXX_

#include <math.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include "Epsilon.h"
#include "LineParametersEstimator.h"

namespace lsqrRecipes {

template< unsigned int dimension >
LineParametersEstimator<dimension>::LineParametersEstimator(double delta) : 
  ParametersEstimator< Point<double, dimension>, double>(2) 
{
  this->deltaSquared = delta*delta;
}
/*****************************************************************************/
/*
 * Estimate the line parameters  [n_0,...,n_k,a_0,...,a_k].
 */
template< unsigned int dimension  >
void LineParametersEstimator<dimension>::estimate(std::vector< Point<double, dimension> *> &data, 
										                              std::vector<double> &parameters)
{
	unsigned int i, j, numParameters;
  double dirNorm;
	
	parameters.clear();
	              //user forgot to initialize the minimal number of required 
                //elements or there are not enough data elements for computation
                //or the points are too close to each other
	if( this->minForEstimate==0 || data.size() < this->minForEstimate ||
      (data[0])->distanceSquared(*(data[1]))< this->deltaSquared )
		return;

  numParameters = 2*dimension;
  parameters.resize(numParameters);
  dirNorm = 0.0;
  for( i=0, j=dimension; i<dimension; i++, j++ ) {
    parameters[i] =  (*data[0])[i] - (*data[1])[i];
    dirNorm += parameters[i]*parameters[i];
    parameters[j] = (*data[0])[i];
  }
  dirNorm = sqrt(dirNorm);
  for( i=0; i<dimension; i++ ) 
    parameters[i]/=dirNorm;
}
/*****************************************************************************/
/*
 * Estimate the line parameters  [n_0,...,n_k,a_0,...,a_k].
 */
template< unsigned int dimension  >
void LineParametersEstimator<dimension>::estimate(std::vector< Point<double, dimension> > &data, 
                                                   std::vector<double> &parameters)
{
	std::vector< Point<double, dimension> *> usedData;
	unsigned int dataSize = static_cast<unsigned int>(data.size());
	for(unsigned int i=0; i<dataSize; i++)
		usedData.push_back(&(data[i]));
	estimate(usedData,parameters);
}
/*****************************************************************************/
/*
 * Estimate the plane parameters  [n_0,...,n_k,a_0,...,a_k].
 */
template< unsigned int dimension  >
void LineParametersEstimator<dimension>::leastSquaresEstimate(std::vector< Point<double, dimension> *> &data, 
                                                              std::vector<double> &parameters)
{
	parameters.clear();	   //not enough data elements for computation
	if(data.size()<this->minForEstimate)
		return;

  unsigned int i, j, k, pointNum = static_cast<unsigned int>(data.size());
  vnl_matrix<double> meanMat(dimension, dimension), covariance(dimension,dimension,0);
  vnl_vector<double> mean(dimension,0);

               //create covariance matrix
  double sqrtN = sqrt((double)pointNum);
  for(i=0; i<pointNum; i++) 
    for(j=0; j<dimension; j++)
      mean[j] += (*data[i])[j];
  mean/= sqrtN;
  for(i=0; i<dimension; i++)
    for(j=i; j<dimension; j++)
       meanMat(i,j) = meanMat(j,i) = mean[i]*mean[j];

                  //upper half
  for(i=0; i<pointNum; i++) 
    for(j=0; j<dimension; j++)
      for(k=j; k<dimension; k++)
        covariance(j,k) += (*data[i])[j] * (*data[i])[k];
                 //copy to lower half
  for(j=0; j<dimension; j++)
    for(k=j+1; k<dimension; k++)
      covariance(k,j) = covariance(j,k);
               //subtract mean matrix
  covariance-=meanMat;

                           //compute eigen-vectors/values of covariance 
  vnl_symmetric_eigensystem<double> eigenSystem(covariance);
  
                //the line direction is the eigenvector corresponding to 
                //the largest eigenvalue
                //I assume ||eigenSystem.V(:,dimension-1)|| = 1
  for(i=0; i<dimension; i++)
    parameters.push_back(eigenSystem.V(i,dimension-1));
  for(i=0; i<dimension; i++)
    parameters.push_back(mean[i]/sqrtN);
}
/*****************************************************************************/
/*
 * Estimate the plane parameters  [n_0,...,n_k,a_0,...,a_k].
 */
template< unsigned int dimension  >
void LineParametersEstimator<dimension>::leastSquaresEstimate(std::vector< Point<double, dimension> > &data, 
														                                  std::vector<double> &parameters)
{
	std::vector< Point<double, dimension> *> usedData;
	unsigned int dataSize = static_cast<unsigned int>(data.size());
	for(unsigned int i=0; i<dataSize; i++)
		usedData.push_back(&(data[i]));
	leastSquaresEstimate(usedData,parameters);
}
/*****************************************************************************/
/*
 * Given the the line parameters  [n_0,...,n_k,a_0,...,a_k] check if
 * ||[data-a]  - dot(data-a,n)n||< delta^2
 *
 * If the parameters vector is too short an exception will
 * be thrown by the vector's [] operator. We don't check vector length explicitly.
 */
template< unsigned int dimension  >
bool LineParametersEstimator<dimension>::agree(std::vector<double> &parameters, 
                                               Point<double, dimension> &data)
{
  double v[dimension], vDotN, distanceSquared;
  unsigned int i, j;

  vDotN = 0.0;
  for(i=0, j=dimension; i<dimension; i++, j++) {
    v[i] = data[i] - parameters[j];
    vDotN+=v[i]*parameters[i];
  }
  distanceSquared = 0.0;
  for(i=0; i<dimension; i++)
    distanceSquared += (v[i] - vDotN*parameters[i])*(v[i] - vDotN*parameters[i]);
	return distanceSquared < this->deltaSquared;
}

} //namespace lsqrRecipes

#endif //_PLANE_PARAMETERS_ESTIMATOR_TXX_
