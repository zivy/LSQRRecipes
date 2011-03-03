#ifndef _ABSOLUTE_ORIENTATION_PARAMETERS_ESTIMATOR_H_
#define _ABSOLUTE_ORIENTATION_PARAMETERS_ESTIMATOR_H_

#include "Copyright.h"

#include "ParametersEstimator.h"
#include "Point3D.h"

/**
 * This class estimates the parameters of a Euclidean (rigid) transformation given corresponding 3D points in
 * two coordinate systems.
 * A Euclidean transformation is represented by a quaterionion and a translation as: 
 *                                              [s,q_x,q_y,q_z,t_x,t_y,t_z]
 *                                              where ||s^2 + q_x^2 + q_y^2 + q_z^2||^2 = 1
 *
 * This implementation is based on the paper "Closed-form solution of absolute orientation using unit quaternions", B.K.P. Horn,
 * Journal of the Optical Society of America, Vol. 4(4), pp 629--642, 1987.
 *
 * @author: Ziv Yaniv (zivy@isis.georgetown.edu)
 *
 */
namespace lsqrRecipes {

class AbsoluteOrientationParametersEstimator : public ParametersEstimator< std::pair<Point3D,Point3D>,double > {
public:
  /**
	 * Object constructor.
	 * @param delta Two points p_1,p_2 correspond via the Euclidean transformation T if ||p_2 - Tp_1||^2 < delta^2
   */
	AbsoluteOrientationParametersEstimator(double delta);

	/**
	 * Compute the Euclidean transformation defined by the given data points such that 
	 * data[i].second = T*data[i].first.
	 * @param data A vector containing three 3D point pairs.
	 * @param This vector is cleared and then filled with the computed parameters, [s,q_x,q_y,q_z,t_x,t_y,t_z].
	 *        If the vector contains less than three point pairs or the points are collinear
	 *        then the resulting parameters vector is empty (size = 0).
	 */
	virtual void estimate(std::vector< std::pair<Point3D,Point3D> *> &data, std::vector<double> &parameters);
	virtual void estimate(std::vector< std::pair<Point3D,Point3D> > &data, std::vector<double> &parameters);
  

	/**
	 * Compute a least squares estimate of the transformation defined by the given point pairs such that
	 * T minimizes sum(||p_2 - Tp_1||^2).
	 * This implementation is of the analytical least squares method presented in 
	 * "Closed-form solution of absolute orientation using unit quaternions", B.K.P. Horn,
   * Journal of the Optical Society of America, Vol. 4(4), pp 629--642, 1987.
	 * @param data The Euclidean transform should minimize the least squares error between all point pairs.
	 * @param parameters This vector is cleared and then filled with the computed parameters, [s,q_x,q_y,q_z,t_x,t_y,t_z].
	 *                   If the vector contains less than three point pairs or all the points are collinear
	 *                   then the resulting parameters vector is empty (size = 0).
	 */
	virtual void leastSquaresEstimate(std::vector< std::pair<Point3D,Point3D> *> &data, std::vector<double> &parameters);
	virtual void leastSquaresEstimate(std::vector< std::pair<Point3D,Point3D> > &data, std::vector<double> &parameters);

	/**
	 * Compute a weighted least squares estimate of the transformation. This is a simple extension to the original 
	 * algorithm described in Horn's paper.
	 * @param data The Euclidean transform should minimize the least squares error between all point pairs.
	 * @param weights Each point pair is weighed according to the given weight. This vector is assumed to have at least
	 *                data.size() entries. Entries are assumed to be non-negative.
	 * @param parameters This vector is cleared and then filled with the computed parameters, [s,q_x,q_y,q_z,t_x,t_y,t_z].
	 *                   If the vector contains less than three point pairs or all the points are collinear
	 *                   then the resulting parameters vector is empty (size = 0).
	 */
	void weightedLeastSquaresEstimate(std::vector< std::pair<Point3D,Point3D> *> &data, std::vector<double> &weights,
		                                std::vector<double> &parameters);


	/**
	 * Return true if the squared distance between mapped point and its matching point
	 * is smaller than 'delta' (see constructor):
	 *   ||(data.second - T(parameters)*data.first)||^2 < delta^2
	 *
	 * @param parameters The parameters (unit quaternion and translation) of the transformation such 
	 *                   that data[i].second = T*data[i].first
	 * @param data Matching points in the two coordinate systems.
	 */
	virtual bool agree(std::vector<double> &parameters, std::pair<Point3D,Point3D> &data);

	/**
	 * Change the distance defining if two points correspond via a given Euclidean transformation.
	 * @param delta The points, p1, p2, correspond if ||p_2 - Tp_1||^2 < delta^2
	 */
	void setDelta(double delta) {this->deltaSquared = delta*delta;}

private:
	double deltaSquared; //given transformation T([s,q_x,q_y,q_z,t_x,t_y,t_z]) two points, p1, p2, correspond if ||p_2 - Tp_1||^2 < delta^2
};

} //namespace lsqrRecipes

#endif //_ABSOLUTE_ORIENTATION_PARAMETERS_ESTIMATOR_H_
