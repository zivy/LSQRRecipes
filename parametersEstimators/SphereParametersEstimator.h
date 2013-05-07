#ifndef _SPHERE_PARAMETERS_ESTIMATOR_H_
#define _SPHERE_PARAMETERS_ESTIMATOR_H_

#include "copyright.h"

#include <vnl/vnl_least_squares_function.h>
#include "ParametersEstimator.h"
#include "Point.h"


/**
 * This class estimates the parameters of a hypersphere (in 2D this is a circle,
 * in 3D this is a sphere...).
 * A sphere is represented as: (1) (p-c)^T(p-c) = sum(p_i-c_i)^2 = r^2
 * where p in R^n is a point on the hypersphere, c in R^n is the sphere's 
 * center, and r its radius.
 * 
 * @author: Ziv Yaniv
 *
 */
namespace lsqrRecipes {

template< unsigned int dimension >
class SphereParametersEstimator : 
  public ParametersEstimator< Point<double, dimension>,double> {
public:

	enum LeastSquaresType {ALGEBRAIC = 0, GEOMETRIC};

	/**
	 * Object constructor.
	 * @param delta A point is on the sphere if its distance from the sphere is 
   *              less than 'delta'.
	 * @param lsType When the leastSquaresEstimate() method is called it computes 
   *               an algebraic or geometric fit. This flag tells it which one.
   */
	SphereParametersEstimator(double delta, LeastSquaresType lsType = GEOMETRIC);

	/**
	 * Compute the hypersphere defined by the given data points.
	 * @param data A vector containing k+1 kD points.
	 * @param parameters This vector is cleared and then filled with the computed 
   *                   parameter values. The parameters of the (hyper)sphere 
   *                   passing through these points [c_0,...,c_k,r].
	 *                   If the vector contains less than k+1 points or the points 
   *                   are linearly dependent (e.g. four coplanar points in the 
   *                   3D case) then the resulting parameters vector is empty 
   *                   (size = 0).
	 */
	virtual void estimate(std::vector< Point<double, dimension> * > &data, 
                        std::vector<double> &parameters);
	virtual void estimate(std::vector< Point<double, dimension> > &data, 
                        std::vector<double> &parameters);

	/**
	 * Compute a least squares estimate of the hypersphere defined by the given 
   * points. This may be either an algebraic or geometric least squares 
   * computation depending on the LeastSquaresType settings. For further 
   * information see geometricLeastSquaresEstimate(), and 
   * algebraicLeastSquaresEstimate().
	 * @param data The sphere should minimize the least squares error to these 
   *             points.
	 * @param parameters This vector is cleared and then filled with the computed 
   *                   parameters. The parameters of the (hyper)sphere passing 
   *                   through these points [c_0,...,c_k,r]. If the vector 
   *                   contains less than k+1 points or the points are linearly 
   *                   dependent (e.g. all points are coplanar in the 3D case) 
   *                   then the resulting parameters vector is empty (size = 0).
	 */
	virtual void leastSquaresEstimate(std::vector< Point<double, dimension> *> &data, 
                                    std::vector<double> &parameters);
	virtual void leastSquaresEstimate(std::vector< Point<double, dimension> > &data, 
                                    std::vector<double> &parameters);

	/**
	 * Return true if the geometric distance between the hypersphere and the 
   * given point is smaller than 'delta' (see constructor).
	 *
	 * @param parameters The sphere parameters [c_0,...,c_k,r].
	 * @param data Check that the geometric distance between this point and the 
   *             sphere is smaller than 'delta'.
   * Note:If the parameters vector is too short the method 
   *      will throw an exception as it tries to access a non existing entry. 
	 */
	virtual bool agree(std::vector<double> &parameters, 
                     Point<double, dimension> &data);

	/**
	 * Change the distance defining if a point is on the hypersphere or not.
	 * @param delta A point is on the sphere if its distance from it is less than
   *             'delta'.
   */	
	void setDelta(double delta) {this->delta = delta;}

	/**
	 * Change the type of least squares solution.
	 * @param lsType When the leastSquaresEstimate() method is called it computes 
   *               an algebraic or geometric fit. This flag tells it which one.
   */
	void setLeastSquaresType(LeastSquaresType lsType) {this->lsType = lsType;}

	/**
	 * Compute a least squares estimate of the hypersphere defined by the given 
   * points.
	 * This implementation is of an algebraic least squares error:
	 * min||Ax-b||, where A = [-2p_0,-2p_1,...,-2p_{k-1},1], 
   *                    x = [c_0,c_1,...,c_{k-1},d], 
   *                    b = [-p_0^2 - p_1^2 - ... -p_{k-1}^2]
	 *
	 * @param data The hypersphere should minimize the algebraic least squares 
   *             error to these points.
	 * @param parameters This vector is cleared and then filled with the computed 
   *                   parameters. The parameters of the (hyper)sphere passing 
   *                   through these points [c_0,...,c_k,r]. If the vector 
   *                   contains less than k+1 points or the points are linearly 
   *                   dependent (e.g. all points are coplanar in the 3D case) 
   *                   then the resulting parameters vector is empty (size = 0).
	 */
	void algebraicLeastSquaresEstimate(std::vector< Point<double, dimension> *> &data, 
                                     std::vector<double> &parameters);
	
  /**
	 * Compute a least squares estimate of the circle defined by the given points.
	 * This implementation is of a geometric least squares error:
	 *    min (sum (r - sqrt((p_0-c_0)^2 + (p_1-c_1)^2 +...+ (p_{k-1}-c_{k-1})^2))^2
	 *
	 * @param data The (hyper)sphere should minimize the geometric least squares 
   *             error to these points.
	 * @param initialParameters This vector contains the initial parameter values 
   *                          for nonlinear estimation.
   * @param finalParameters This vector is cleared and then filled with the 
   *                        computed sphere parameters [centerX,centerY,centerZ,r].
   *
   * NOTE: As this method encapsulates nonlinear optimization (Levenberg-Marquardt) 
   *       it may not converge due to the specific choice of LM parameter settings.
   *       Possibly modify these values to be appropriate for user's specific input.
	 */
	void geometricLeastSquaresEstimate(std::vector< Point<double, dimension> *> &data, 
																		 std::vector<double> &initialParameters, 
                                     std::vector<double> &finalParameters);

  /**
   * Get the distances between the given hypersphere and the given points.
   * These are pushed at the end of the 'distances' vector.
   *
   * @param parameters The (hyper)sphere parameters [c,r].
   * @param data The points whose distance from the (hyper)sphere we want.   
   * @param distances The computed distances. If the vector is empty then the 
   *                  point and distance indexes match,
   *                  otherwise there is an offset due to the original number
   *                  of entries previously found in the vector.
   * @param min Minimal distance.
   * @param max Maximal distance.
   * @param mean Mean distance.
   */
  static void getDistanceStatistics( std::vector<double> &parameters, 
                                     std::vector< Point<double, dimension> > &data, 
                                     std::vector<double> &distances,
                                     double &min, double &max, double &mean );

private:
  enum {CIRCLE = 2, SPHERE=3};
	double delta; 
            //algebraic or geometric least squares
	LeastSquaresType lsType; 

	static const double SPHERE_EPS;

  inline void estimate2D(std::vector< Point<double, dimension> * > &data, 
                         std::vector<double> &parameters);
  inline void estimate3D(std::vector< Point<double, dimension> * > &data, 
                         std::vector<double> &parameters);
  inline void estimateND(std::vector< Point<double, dimension> * > &data, 
                         std::vector<double> &parameters);

  class SumSquaresSpherePointsDistanceFunction : public vnl_least_squares_function {
		public:
			SumSquaresSpherePointsDistanceFunction(std::vector<Point<double, dimension> *> *data);
			void setPoints(std::vector<Point<double, dimension> *> *data);
      virtual void f( vnl_vector<double> const &x, vnl_vector<double> &fx );
      virtual void gradf( vnl_vector<double> const& x, vnl_matrix<double>& jacobian );
		private:
			std::vector<Point<double, dimension> *> *data;
	};
	
};

} //namespace lsqrRecipes

#include "SphereParametersEstimator.hxx" //the implementation is in this file

#endif //_SPHERE_PARAMETERS_ESTIMATOR_H_
