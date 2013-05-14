#include <math.h>
#include "Line2DParametersEstimator.h"

namespace lsqrRecipes {

Line2DParametersEstimator::Line2DParametersEstimator(double delta) : ParametersEstimator<Point2D,double>(2) {this->deltaSquared = delta*delta;}
/*****************************************************************************/
/*
 * Compute the line parameters  [n_x,n_y,a_x,a_y]
 */
void Line2DParametersEstimator::estimate(std::vector<Point2D *> &data, 
																	std::vector<double> &parameters)
{
	parameters.clear();
	              //not enough data elements for computation
	if(data.size()<this->minForEstimate)
		return;

	double nx = (*data[1])[1] - (*data[0])[1];
	double ny = (*data[0])[0] - (*data[1])[0];
	double normSquared = nx*nx + ny*ny;
	              //the two points are too close to each other
	if(normSquared<this->deltaSquared)
		return;
	double norm = sqrt(nx*nx + ny*ny);
	

	parameters.push_back(nx/norm);
	parameters.push_back(ny/norm);
	parameters.push_back((*data[0])[0]);
	parameters.push_back((*data[0])[1]);		
}
/*****************************************************************************/
/*
 * Compute the line parameters  [n_x,n_y,a_x,a_y]
 */
void Line2DParametersEstimator::estimate(std::vector<Point2D> &data, 
																	std::vector<double> &parameters)
{
	std::vector<Point2D *> usedData;
	int dataSize = data.size();
	for(int i=0; i<dataSize; i++)
		usedData.push_back(&(data[i]));
	estimate(usedData,parameters);
}
/*****************************************************************************/
/*
 * Compute the line parameters  [n_x,n_y,a_x,a_y]
 */
void Line2DParametersEstimator::leastSquaresEstimate(std::vector<Point2D *> &data, 
																							std::vector<double> &parameters)
{
	double meanX, meanY, nx, ny, norm;
	double covMat11, covMat12, covMat21, covMat22; // The entries of the symmetric covariance matrix
	int i, dataSize = data.size();

	parameters.clear();
	if(data.size()<this->minForEstimate)
		return;

	meanX = meanY = 0.0;
	covMat11 = covMat12 = covMat21 = covMat22 = 0;
	for(i=0; i<dataSize; i++) {
		meanX +=(*data[i])[0];
		meanY +=(*data[i])[1];

		covMat11	+=(*data[i])[0] * (*data[i])[0];
		covMat12	+=(*data[i])[0] * (*data[i])[1];
		covMat22	+=(*data[i])[1] * (*data[i])[1];
	}

	meanX/=dataSize;
	meanY/=dataSize;

	covMat11 -= dataSize*meanX*meanX;
  covMat12 -= dataSize*meanX*meanY;
	covMat22 -= dataSize*meanY*meanY;
	covMat21 = covMat12;

	if(covMat11<1e-12) {
		nx = 1.0;
	  ny = 0.0;
		if(covMat22<1e-12)  //both diagonal entries of the covariance matrix are small (variance in x and y directions)
			return;           //which means all the points are the "same" point
	}
	else {	    //lambda1 is the largest eigenvalue of the covariance matrix 
	           //and is used to compute the eigenvector corresponding to the smallest
	           //eigenvalue, which isn't computed explicitly.
		double lambda1 = (covMat11 + covMat22 + sqrt((covMat11-covMat22)*(covMat11-covMat22) + 4*covMat12*covMat12)) / 2.0;
		nx = -covMat12;
		ny = lambda1 - covMat22;
		norm = sqrt(nx*nx + ny*ny);
		nx/=norm;
		ny/=norm;
	}
	parameters.push_back(nx);
	parameters.push_back(ny);
	parameters.push_back(meanX);
	parameters.push_back(meanY);
}
/*****************************************************************************/
/*
 * Compute the line parameters  [n_x,n_y,a_x,a_y]
 */
void Line2DParametersEstimator::leastSquaresEstimate(std::vector<Point2D> &data, 
																							std::vector<double> &parameters)
{
	std::vector<Point2D *> usedData;
	int dataSize = data.size();
	for(int i=0; i<dataSize; i++)
		usedData.push_back(&(data[i]));
	leastSquaresEstimate(usedData,parameters);
}
/*****************************************************************************/
/*
 * Given the line parameters  [n_x,n_y,a_x,a_y] check if
 * [n_x, n_y] dot [data.x-a_x, data.y-a_y] < delta
 */
bool Line2DParametersEstimator::agree(std::vector<double> &parameters, Point2D &data)
{
	double signedDistance = parameters[0]*(data[0]-parameters[2]) + parameters[1]*(data[1]-parameters[3]); 
	return ((signedDistance*signedDistance) < this->deltaSquared);
}

} //namespace lsqrRecipes