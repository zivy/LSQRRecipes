#ifndef _PLANE_PARAMETERS_ESTIMATOR_H_
#define _PLANE_PARAMETERS_ESTIMATOR_H_

#include "copyright.h"

#include "ParametersEstimator.h"
#include "Point.h"

/**
 * This class estimates the parameters of a (hyper)plane.
 * This is a generalization of the 2D line case. The reason there is also a class for the 2D line case 
 * is that in the 2D case we do not need a linear algebra package for computing 
 * eigen-vectors. 
 *
 * 
 * A (hyper)plane is represented as: (*) dot(n,p-a)=0 
 *                                    where n is the (hyper)plane normal (|n| = 1) and 'a' is a 
 *                                    point on the (hyper)plane. 
 * All points 'p' which satisfy equation (*) are on the (hyper)plane.
 * 
 * @author: Ziv Yaniv (zivy@isis.georgetown.edu)
 *
 */
namespace lsqrRecipes {

template< unsigned int dimension >
class PlaneParametersEstimator : public ParametersEstimator< Point<double, dimension>,double> {
public:
  /**
	 * Object constructor.
	 * @param delta A point is on the (hyper)plane if its distance from the hyper(plane) is less than 'delta'.
   */
	PlaneParametersEstimator(double delta);

	/**
	 * Compute the (hyper)plane defined by the given data points.
	 * @param data A vector containing k kD points.
	 * @param This vector is cleared and then filled with the computed parameters.
	 *        The parameters of the plane passing through these points [n_0,...,n_k,a_0,...,a_k]
	 *        where ||(n_0,...,nk)|| = 1.
	 *        If the vector contains less than k points or the points are linearly dependent 
	 *        then the resulting parameters vector is empty (size = 0).
	 */
	virtual void estimate(std::vector< Point<double, dimension> *> &data, std::vector<double> &parameters);
	virtual void estimate(std::vector< Point<double, dimension> > &data, std::vector<double> &parameters);

	/**
	 * Compute a least squares estimate of the line defined by the given points.
	 * This implementation is of an orthogonal least squares error.
	 *
	 * @param data The line should minimize the least squares error to these points.
	 * @param parameters This vector is cleared and then filled with the computed parameters.
	 *                   Fill this vector with the computed line parameters [n_x,n_y,a_x,a_y]
	 *                   where ||(n_x,ny)|| = 1.
	 *                   If the vector contains less than two points or all the points are coincident 
	 *                   then the resulting parameters vector is empty (size = 0).
	 */
	virtual void leastSquaresEstimate(std::vector< Point<double, dimension> *> &data, std::vector<double> &parameters);
	virtual void leastSquaresEstimate(std::vector< Point<double, dimension> > &data, std::vector<double> &parameters);

	/**
	 * Return true if the distance between the line defined by the parameters and the
	 * given point is smaller than 'delta' (see constructor).
	 * @param parameters The line parameters [n_x,n_y,a_x,a_y].
	 * @param data Check that the distance between this point and the line is smaller than 'delta'.
   *
   * Note:If the parameters vector is too short the method will throw an 
   *      exception as it tries to access a non existing entry. 
   *
	 */
	virtual bool agree(std::vector<double> &parameters, Point<double, dimension> &data);

	/**
	 * Change the distance defining if a point is on the line or not.
	 * @param delta A point is on the line if its distance from it is less than 'delta'.
	 */
	void setDelta(double delta) {this->deltaSquared = delta*delta;}

private:
	double deltaSquared; //given line L and point P, if dist(L,P)^2 < delta^2 then the point is on the line
};

} //namespace lsqrRecipes

#include "PlaneParametersEstimator.txx" //the implementation is in this file

#endif //_PLANE_PARAMETERS_ESTIMATOR_H_
