#ifndef _LINE_PARAMETERS_ESTIMATOR_H_
#define _LINE_PARAMETERS_ESTIMATOR_H_

#include "copyright.h"

#include "ParametersEstimator.h"
#include "Point.h"

/**
 * This class estimates the parameters of a line in R^d.
 * 
 * A line is represented as the set of points: (*) a+tn 
 *                                    where 'a' is a point on the line,
 *                                    't' is a real number, and 
 *                                    n is a vector which is parallel to 
 *                                    the line (|n| = 1).
 *
 * A point 'p' is on the line if it satisfies the following equation:
 *                  ||(p-a) - n^T(p-a)n|| = 0
 *
 * Note: A hyperplane (see PlaneParametersEstimator) in R^2 is a line. 
 *       Hyperplane representation: [n,a] 'n' is the line normal and 'a' is a 
 *                                  point on the line and.
 *       Line representation: [n,a] 'n' is the line's direction (perpendicular 
 *                            to the hyperplane's 'n') and 'a' is a point on the 
 *                            line.
 *        
 * 
 * @author: Ziv Yaniv 
 *
 */
namespace lsqrRecipes {

template< unsigned int dimension >
class LineParametersEstimator : public ParametersEstimator< Point<double, dimension>,double> {
public:
  /**
	 * Object constructor.
	 * @param delta A point is on the line if its distance from the line is less 
   *              than 'delta'.
   */
	LineParametersEstimator(double delta);

	/**
	 * Compute the line defined by the given data points.
	 * @param data A vector containing two kD points.
	 * @param This vector is cleared and then filled with the computed parameters.
	 *        The parameters of the line passing through these points 
   *        [n_0,...,n_k,a_0,...,a_k] where ||(n_0,...,nk)|| = 1.
	 *        If the vector contains less than two points or the points are 
   *        too close (distance between them is less than 'delta')
	 *        then the resulting parameters vector is empty (size = 0).
	 */
	virtual void estimate(std::vector< Point<double, dimension> *> &data, 
                        std::vector<double> &parameters);
	virtual void estimate(std::vector< Point<double, dimension> > &data, 
                        std::vector<double> &parameters);

	/**
	 * Compute a least squares estimate of the line defined by the given points.
	 * This implementation is of an orthogonal least squares error.
	 *
	 * @param data The line should minimize the least squares error to these points.
	 * @param parameters This vector is cleared and then filled with the computed parameters.
	 *                   Fill this vector with the computed line parameters [n,a]
	 *                   where ||n|| = 1.
	 *                   If the vector contains less than two points or all the points are coincident 
	 *                   then the resulting parameters vector is empty (size = 0).
	 */
	virtual void leastSquaresEstimate(std::vector< Point<double, dimension> *> &data, 
                                    std::vector<double> &parameters);
	virtual void leastSquaresEstimate(std::vector< Point<double, dimension> > &data, 
                                    std::vector<double> &parameters);

	/**
	 * Return true if the distance between the line defined by the parameters and the
	 * given point is smaller than 'delta' (see constructor).
	 * @param parameters The line parameters [n,a].
	 * @param data Check that the distance between this point and the line is smaller than 'delta'.
     *
     * Note:If the parameters vector is too short the method will throw an 
     *      exception as it tries to access a non existing entry. 
     *
	 */
	virtual bool agree(std::vector<double> &parameters, Point<double, dimension> &data);

	/**
	 * Change the distance defining if a point is on the line or not.
	 * @param delta A point is on the line if its distance from it is less 
   *              than 'delta'.
	 */
	void setDelta(double delta) {this->deltaSquared = delta*delta;}

private:
    //given line L and point P, if dist(L,P)^2 < delta^2 then the point 
    //is on the line
	double deltaSquared; 
};

} //namespace lsqrRecipes

#include "LineParametersEstimator.hxx" //the implementation is in this file

#endif //_LINE_PARAMETERS_ESTIMATOR_H_
