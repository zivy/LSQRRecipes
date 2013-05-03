#ifndef _PLANE_PARAMETERS_ESTIMATOR_H_
#define _PLANE_PARAMETERS_ESTIMATOR_H_

#include "copyright.h"

#include "ParametersEstimator.h"
#include "Point.h"

/**
 * This class estimates the parameters of a hyperplane (in 2D this is a line, 
 * in 3D a plane...).
 * 
 * A hyperplane is represented as: (*) dot(n,p-a)=0 
 *                                    where n is the hyperplane normal (|n| = 1) and 'a' is a 
 *                                    point on the hyperplane. 
 * All points 'p' which satisfy equation (*) are on the hyperplane.
 * 
 * @author: Ziv Yaniv
 *
 */
namespace lsqrRecipes {

template< unsigned int dimension >
class PlaneParametersEstimator : public ParametersEstimator< Point<double, dimension>,double> {
public:
  /**
	 * Object constructor.
	 * @param delta A point is on the hyperplane if its distance from the 
   *              hyperplane is less than 'delta'.
   */
	PlaneParametersEstimator(double delta);

	/**
	 * Compute the hyperplane defined by the given data points.
	 * @param data A vector containing k kD points.
	 * @param This vector is cleared and then filled with the computed parameters.
	 *        The parameters of the plane passing through these points 
   *        [n_0,...,n_k,a_0,...,a_k] where ||(n_0,...,n_k)|| = 1.
	 *        If the vector contains less than k points or the points are linearly 
   *        dependent then the resulting parameters vector is empty (size = 0).
	 */
  virtual void estimate(std::vector< Point<double, dimension> *> &data, 
                        std::vector<double> &parameters);
  virtual void estimate(std::vector< Point<double, dimension> > &data, 
                        std::vector<double> &parameters);

	/**
	 * Compute a least squares estimate of the hyperplane defined by the given points.
	 * This implementation is of an orthogonal least squares error.
	 *
	 * @param data The hyperplane should minimize the least squares error to these 
   *             points.
	 * @param parameters This vector is cleared and then filled with the computed 
   *                   parameters. Fill this vector with the computed hyperplane 
   *                   parameters [n_0,...,n_k,a_0,...,a_k] where ||(n_0,...,n_k)|| = 1.
	 *                   If the vector contains less than two points or all the 
   *                   points are coincident then the resulting parameters 
   *                   vector is empty (size = 0).
	 */
  virtual void leastSquaresEstimate(std::vector< Point<double, dimension> *> &data, 
                                    std::vector<double> &parameters);
  virtual void leastSquaresEstimate(std::vector< Point<double, dimension> > &data, 
                                    std::vector<double> &parameters);

	/**
	 * Return true if the distance between the hyperplane defined by the parameters and 
   * the given point is smaller than 'delta' (see constructor).
	 * @param parameters The hyperplane parameters [n_0,...,n_k,a_0,...,a_k].
	 * @param data Check that the distance between this point and the hyperplane is 
   *             less than 'delta'.
   *
   * Note:If the parameters vector is too short the method will throw an 
   *      exception as it tries to access a non existing entry. 
   *
	 */
	virtual bool agree(std::vector<double> &parameters, Point<double, dimension> &data);

	/**
	 * Change the distance defining if a point is on the hyperplane or not.
	 * @param delta A point is on the hyperplane if its distance from it is less 
   *              than 'delta'.
	 */
	void setDelta(double delta) {this->deltaSquared = delta*delta;}

private:
    //given hyperplane H(a,n) and point P, if dist(H,P)^2 < delta^2 then the 
    //point is on the hyperplane
	double deltaSquared; 
};

} //namespace lsqrRecipes

#include "PlaneParametersEstimator.txx" //the implementation is in this file

#endif //_PLANE_PARAMETERS_ESTIMATOR_H_
