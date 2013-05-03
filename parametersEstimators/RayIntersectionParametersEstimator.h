#ifndef _RAY_INTERSECTION_PARAMETERS_ESTIMATOR_H_
#define _RAY_INTERSECTION_PARAMETERS_ESTIMATOR_H_

#include "copyright.h"

#include "ParametersEstimator.h"
#include "Point3D.h"
#include "Ray3D.h"

/**
 * This class estimates the intersection point of a given set of three 
 * dimensional rays.
 *
 * NOTE: It is assumed that the ray direction has unit norm (||n||=1). This 
 *       assumption is not enforced by the Ray3D class. The user is expected to 
 *       conform (i.e. normalize the ray direction).
 *
 * @author: Ziv Yaniv
 *
 */

namespace lsqrRecipes {

class RayIntersectionParametersEstimator : public ParametersEstimator<Ray3D,double> {
public:
  /**
	 * Object constructor.
	 * @param delta A ray is coincident with an intersection point if the 
   *              point-ray distance is less than 'delta'.
   * @param minimalAngularDeviation Two rays are considered parallel if the 
   *                                angular deviation between them is less than 
   *                                this value [radians]. Default is one degree.
   */
	RayIntersectionParametersEstimator(double delta, 
                                     double minimalAngularDeviation = 0.017453292519943295769236907684886);

	/**
	 * Compute the intersection point defined by the given rays. This is the 
   * unique point that minimizes the sum of distances to both rays.
	 * @param data A vector containing two rays.
	 * @param parameters This vector is cleared and then filled with the computed 
   *                   intersection point parameters [x,y,z]. If the vector 
   *                   contains less than two rays or they (a) do not intersect 
   *                   or (b) have an infinite number of intersection points 
   *                   then the resulting parameters vector is empty (size = 0).
	 */
	virtual void estimate(std::vector<Ray3D *> &data, 
                        std::vector<double> &parameters);
	virtual void estimate(std::vector<Ray3D> &data, 
                        std::vector<double> &parameters);

	/**
	 * Compute a least squares estimate of the intersection point defined by the 
   * given rays.
	 * NOTE: In practice the rays are treated as lines so the point we get is the 
   *       one that minimizes its distance from all lines {t in (-inf,inf)} and 
   *       not rays {t in [0,inf)}. Use this method only when line intersection 
   *       and ray intersection are equivalent. In general they are not 
   *       equivalent, the following two rays don't intersect but the lines do:
   *
   *                                      *   *------------>
   *                                      | 
   *                                      |
   *                                      \/
   *
   *       When this method is used as part of the RANSAC based estimation we 
   *       are sure that the ray parameter t is in [0,inf) as the methods 
   *       estimate() and agree() enforce this constraint, which in turn means 
   *       that treating these rays as lines will result in the correct answer.
	 * @param data The point should minimize the least squares error (distance) to
   *             these rays.
	 * @param parameters This vector is cleared and then filled with the computed 
   *                   parameters point parameters [x,y,z]. If the vector 
   *                   contains less than two rays or they (a) do not intersect 
   *                   or (b) have an infinite number of intersection points 
   *                   then the resulting parameters vector is empty (size = 0).
	 */
	virtual void leastSquaresEstimate(std::vector<Ray3D *> &data, 
                                    std::vector<double> &parameters);
	virtual void leastSquaresEstimate(std::vector<Ray3D> &data, 
                                    std::vector<double> &parameters);

	/**
	 * Return true if the distance between the point defined by the parameters and
   * the given ray is smaller than 'delta' (see constructor).
	 * @param parameters The point parameters [x,y,z].
	 * @param data Check that the distance between this ray and the point is 
   *             smaller than 'delta'.
   *
   * Note:If the parameters vector contains less than three entries the method 
   *      will throw an exception as it tries to access a non existing entry. 
	 */
	virtual bool agree(std::vector<double> &parameters, Ray3D &data);

	/**
	 * Change the distance defining if a ray and a given point are considered 
   * coincident.
	 * @param delta A ray is coincident with an intersection point if the 
   *              point-ray distance is less than 'delta'.
	 */
	void setDelta(double delta) {this->deltaSquared = delta*delta;}

private:
     //given ray R and point P, if dist(R,P)^2 < delta^2 then the 
     //point is on the ray
	double deltaSquared;
     //crossEps = sin(minimalAngularDeviation) two rays with normals n1,n2 are 
     //parallel if cross(n1,n2)<crossEps, assuming ||n1||==||n2|| = 1
  double crossEps; 
};

} //namespace lsqrRecipes

#endif //_RAY_INTERSECTION_PARAMETERS_ESTIMATOR_H_
