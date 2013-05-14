#ifndef _LINE2D_PARAMETERS_ESTIMATOR_H_
#define _LINE2D_PARAMETERS_ESTIMATOR_H_

#include "copyright.h"

#include "ParametersEstimator.h"
#include "Point2D.h"

/**
 * This class estimates the parameters of 2D lines.
 * It uses exactly the same approach as that used for hyperplane estimation. The
 * only reason to use this class instead of the PlaneParametersEstimator with
 * dimension 2 is that you can easily rip it out of the Least Squares Recipes 
 * library and incorporate it into your codebase (it does not require an
 * external library to compute the eigenvectors as required by the more generic
 * hyperplane estimation class).
 *
 * Notes: The parameterization is the same as that of the hyperplane [n,a] where
 *        'n' is the normal to the line (unique only in 2D), and 'a' is a point
 *        on the line. This is different from the parameterization used for the
 *        kD line where 'n' is the direction parallel to the line.
 *
 *        The unit test for this class is incorporated into the unit test of the
 *        more generic LineParametersEstimator. 
 *
 * Author: Ziv Yaniv
 */
 
namespace lsqrRecipes {

class Line2DParametersEstimator : public ParametersEstimator<Point2D,double> {
public:
  /**
	 * Object constructor.
	 * @param delta A point is on the line if its distance from the line is less 
   *              than 'delta'.
   */
	Line2DParametersEstimator(double delta);

	/**
	 * Compute the line defined by the given data points.
	 * @param data A vector containing two 2D points.
	 * @param This vector is cleared and then filled with the computed parameters.
	 *        The parameters of the line passing through these points [n_x,n_y,a_x,a_y]
	 *        where ||(n_x,ny)|| = 1.
	 *        If the vector contains less than two points or the points are coincident 
	 *        (closer than the user specified delta) then the resulting parameters 
   *        vector is empty (size = 0).
	 */
	virtual void estimate(std::vector<Point2D *> &data, std::vector<double> &parameters);
	virtual void estimate(std::vector<Point2D> &data, std::vector<double> &parameters);

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
	virtual void leastSquaresEstimate(std::vector<Point2D *> &data, std::vector<double> &parameters);
	virtual void leastSquaresEstimate(std::vector<Point2D> &data, std::vector<double> &parameters);

	/**
	 * Return true if the distance between the line defined by the parameters and the
	 * given point is smaller than 'delta' (see constructor).
	 * @param parameters The line parameters [n_x,n_y,a_x,a_y].
	 * @param data Check that the distance between this point and the line is smaller than 'delta'.
	 */
	virtual bool agree(std::vector<double> &parameters, Point2D &data);

	/**
	 * Change the distance defining if a point is on the line or not.
	 * @param delta A point is on the line if its distance from it is less than 'delta'.
	 */
	void setDelta(double delta) {this->deltaSquared = delta*delta;}

private:
    //given line L and point P, if dist(L,P)^2 < delta^2 then the point is on the line
	double deltaSquared; 
};

} //namespace lsqrRecipes

#endif //_LINE2D_PARAMETERS_ESTIMATOR_H_