#include<vnl/algo/vnl_matrix_inverse.h>
#include "Epsilon.h"
#include "PivotCalibrationParametersEstimator.h"

namespace lsqrRecipes {

PivotCalibrationEstimator::PivotCalibrationEstimator(double delta) : ParametersEstimator< Frame, double >(3) {this->delta = delta;}
/*****************************************************************************/
void PivotCalibrationEstimator::estimate(std::vector< Frame *> &data, 
                												 std::vector<double> &parameters)
{
	parameters.clear();
	              //not enough data elements for computation
	if(data.size()<this->minForEstimate)
		return;

  double minusIValues[] = {-1.0, 0.0, 0.0, 0.0,-1.0, 0.0, 0.0, 0.0, -1.0};
  vnl_matrix<double> A(9,6), R(3,3), minusI(minusIValues,3,3);
  vnl_vector<double> x(6), b(9), t(3);

  
  data[0]->getRotationMatrix(R);
  A.update(R,0,0);
  A.update(minusI,0,3);
  data[0]->getTranslation(t);
  b.update(-t,0);

  data[1]->getRotationMatrix(R);
  A.update(R,3,0);
  A.update(minusI,3,3);
  data[1]->getTranslation(t);
  b.update(-t,3);

  data[2]->getRotationMatrix(R);
  A.update(R,6,0);
  A.update(minusI,6,3);
  data[2]->getTranslation(t);
  b.update(-t,6);

  vnl_matrix_inverse<double> Ainv(A);
                     //explicitly zero out small singular values
                     //this is ugly as it exposes that the inverse is computed via SVD
  Ainv.zero_out_absolute(EPS);

  if(Ainv.rank()<6) //matrix is not invertible
    return;
  x = Ainv * b;

  for(unsigned int i=0; i<6; i++)
    parameters.push_back(x[i]);
}
/*****************************************************************************/
void PivotCalibrationEstimator::estimate(std::vector< Frame > &data, 
																	       std::vector<double> &parameters)
{
	std::vector< Frame *> usedData;
	unsigned int dataSize = static_cast<unsigned int>(data.size());
	for(unsigned int i=0; i<dataSize; i++)
		usedData.push_back(&(data[i]));
	estimate(usedData, parameters);
}
/*****************************************************************************/
void PivotCalibrationEstimator::leastSquaresEstimate(std::vector< Frame *> &data, 
																							       std::vector<double> &parameters)
{
	parameters.clear();
	              //not enough data elements for computation
	if(data.size()<this->minForEstimate)
		return;

  unsigned int n = data.size();
  double minusIValues[] = {-1.0, 0.0, 0.0, 0.0,-1.0, 0.0, 0.0, 0.0, -1.0};
  vnl_matrix<double> A(3*n,6), R(3,3), minusI(minusIValues,3,3);
  vnl_vector<double> x(6), b(3*n), t(3);

  
  for(unsigned int i=0; i<n; i++) {
    data[i]->getRotationMatrix(R);
    A.update(R,3*i,0);
    A.update(minusI,3*i,3);
    data[i]->getTranslation(t);
    b.update(-t,3*i);
  }

  vnl_matrix_inverse<double> Ainv(A);
                     //explicitly zero out small singular values
                     //this is ugly as it exposes that the inverse is computed via SVD
  Ainv.zero_out_absolute(EPS);

  if(Ainv.rank()<6) //we don't have a unique solution to the normal equations
    return;
  x = Ainv * b;
  for(unsigned int i=0; i<6; i++)
    parameters.push_back(x[i]);

}
/*****************************************************************************/
void PivotCalibrationEstimator::leastSquaresEstimate(std::vector< Frame > &data, 
																							       std::vector<double> &parameters)
{
	std::vector< Frame *> usedData;
	unsigned int dataSize = static_cast<unsigned int>(data.size());
	for(unsigned int i=0; i<dataSize; i++)
		usedData.push_back(&(data[i]));
	leastSquaresEstimate(usedData, parameters);
}
/*****************************************************************************/
bool PivotCalibrationEstimator::agree(std::vector<double> &parameters, Frame &data)
{
  Point3D tDRF, tW;

  tDRF[0] = parameters[0];
  tDRF[1] = parameters[1];
  tDRF[2] = parameters[2];
  tW[0] = parameters[3];
  tW[1] = parameters[4];
  tW[2] = parameters[5];

         //apply the transformation to the point, result is written in place
  data.apply(tDRF);
            //distance between the two points is less than the threshold
  return (tDRF - tW).l2Norm()<this->delta;
}

} //namespace lsqrRecipes
