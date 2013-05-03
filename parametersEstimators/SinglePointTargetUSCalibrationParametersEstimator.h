#ifndef _SINGLE_POINT_TARGET_US_CALIBRATION_PARAMETERS_ESTIMATOR_H_
#define _SINGLE_POINT_TARGET_US_CALIBRATION_PARAMETERS_ESTIMATOR_H_

#include "copyright.h"

#include <vnl/vnl_least_squares_function.h>
#include "ParametersEstimator.h"
#include "Point3D.h"
#include "Point2D.h"
#include "Frame.h"


/**
 * These two classes estimate the transformation mapping points from the US image
 * coordinate system to the coordinate system of a tracked reference frame that
 * is fixed relative to the US probe. The method is based on a single point 
 * target. This can be either a cross wire phantom or a phantomless approach, 
 * using a tracked pointer. The two approaches are very similar but are not 
 * exactly the same. When using a cross wire phantom the location of the point 
 * target in the phantom's coordinate system is estimated as part
 * of the problem formulation. When using a calibrated pointer, the location of
 * the pointer tip relative to the pointer's dynamic reference frame is known.
 * 
 * For additional information with regard to US 
 * calibration see: 
 *
 * "A review of calibration techniques for freehand 3-D ultrasound systems", 
 * L. Mercier, T. Lango, F. Lindseth, L.D. Collins, Ultrasound in Med. & Biol., 
 * vol. 31(2), pp. 143-165, 2005.
 * 
 * @author: Ziv Yaniv
 *
 */

namespace lsqrRecipes {

/**
 * The data elements used by the phantom based estimator are pairs that 
 * include:
 * (a). T2 - the transformation mapping points from the US dynamic reference 
 *           frame to the fixed, tracker, coordinate system T.
 * (b). q - the target point's coordinates in the image coordinate system 
 *          [u,v,0].
 */
struct SingleUnknownPointTargetUSCalibrationParametersEstimatorDataType {
    Frame T2;  
    Point2D q;
};


class SingleUnknownPointTargetUSCalibrationParametersEstimator : 
  public ParametersEstimator<SingleUnknownPointTargetUSCalibrationParametersEstimatorDataType,double> {
public:
  
  enum LeastSquaresType {ANALYTIC = 0, ITERATIVE};

  typedef SingleUnknownPointTargetUSCalibrationParametersEstimatorDataType DataType;

  /**
	 * Object constructor.
	 * @param delta A data element (pair of US-DRF transformation and 2D point location) 
   *              agrees with specific calibration parameter values if
   *              the 3D distance between the mapped point and its known location
   *              is less than 'delta'.
	 * @param lsType When the leastSquaresEstimate() method is called it computes 
   *               an analytic or iterative fit. This flag tells it which one.
   */
	SingleUnknownPointTargetUSCalibrationParametersEstimator(double delta, 
                                                           LeastSquaresType lsType = ITERATIVE);

	/**
	 * Compute the affine US calibration transformation defined by the given data 
   * elements. 
   *
   * @param data A vector containing four data elements.
	 * @param parameters This vector is cleared and then filled with the computed 
   *                   transformation parameters 
   * [t_1x, t_1y, t_1z, t_3x, t_3y, t_3z, omega_z, omega_y, omega_x, m_x, m_y, m_x*R3(:,1), m_y*R3(:,2), R3(:,3)].
   *
   *                   Angles are in radians.
   *
   *                   The reason for providing the composition of rotation and 
   *                   scale, last nine entries, is computational efficiency of
   *                   the agree method. 
   *
   *                   If the vector contains less than four data elements or 
   *                   they do not provide a set of linearly independent equation
   *                   systems the resulting parameters vector is empty 
   *                   (size == 0).
	 */	
  virtual void estimate(std::vector<DataType *> &data, 
                        std::vector<double> &parameters);
	virtual void estimate(std::vector<DataType> &data, 
                        std::vector<double> &parameters);


	/**
	 * Compute the least squares affine US calibration transformation defined by 
   * the given data elements. This may be either an analytic or iterative least 
   * squares computation depending on the LeastSquaresType settings. For further 
   * information see iterativeLeastSquaresEstimate(), and 
   * analyticLeastSquaresEstimate().
   *
   * @param data A vector containing at least four data elements.
	 * @param parameters This vector is cleared and then filled with the computed 
   *                   transformation parameters 
   * [t_1x, t_1y, t_1z, t_3x, t_3y, t_3z, omega_z, omega_y, omega_x, m_x, m_y, m_x*R3(:,1), m_y*R3(:,2), R3(:,3)].
   *
   *                   Angles are in radians.
   *
   *                   The reason for providing the composition of rotation and 
   *                   scale, last nine entries, is computational efficiency of
   *                   the agree method. 
   *
   *                   If the vector contains less than four data elements or 
   *                   they do not provide a set of linearly independent equation
   *                   systems the resulting parameters vector is empty 
   *                   (size == 0).
	 */	
	virtual void leastSquaresEstimate(std::vector<DataType *> &data, 
                                    std::vector<double> &parameters);
	virtual void leastSquaresEstimate(std::vector<DataType> &data, 
                                    std::vector<double> &parameters);

  /**
   * A given data element is consistent, agrees, with the parameter values if 
   * the distance between the mapped 3D point and its known location is less 
   * than a user set threshold (see setDelta()).
   *
	 * @param parameters This vector contains the US calibration transformation 
   *                   parameters and the estimtated 3D point location in the
   *                   world/tracker coordinate system (t_1).
   * [t_1x, t_1y, t_1z, t_3x, t_3y, t_3z, omega_z, omega_y, omega_x, m_x, m_y, m_x*R3(:,1), m_y*R3(:,2), R3(:,3)].
   *
   * Angles are in radians.
   *
   * @param data A data element comprised of the US-DRF transformation and the
   *             2D point location in the image coordinate system.
   *
   * Note:If the parameters vector is too short the method will throw an 
   *      exception as it tries to access a non existing entry. 
   *
   */
	virtual bool agree(std::vector<double> &parameters, DataType &data);

  /**
	 * Change the type of least squares solution.
	 * @param lsType When the leastSquaresEstimate() method is called it computes 
   *               an analytic or iterative fit. This flag tells it which one.
   */
	void setLeastSquaresType(LeastSquaresType lsType) {this->lsType = lsType;}

  /**
   * Set the threshold that defines if a data element agrees with a specific 
   * model estimate.
	 * @param delta A data element (pair of US-DRF transformation and 2D point location) 
   *              agrees with specific calibration parameter values if
   *              the 3D distance between the mapped point and its known location
   *              is less than 'delta'.
   */
	void setDelta(double delta) {this->deltaSquared = delta*delta;}

	/**
	 * Compute a least squares estimate of the US calibration and target point
   * location in the world/tracker coordinate system, as defined by the given 
   * data elements (pairs of US-DRF transformation and 2D point location).
	 * This implementation is of an analytic least squares error:
	 * min||Ax-b||, where A = [u_i*R_2i, v_i*R_2i, R_2i -I], 
   *                    x = [m_x*R_3(:,1), m_y*R_3(:,2), t_3, t_1], 
   *                    b = [-t_2i]
	 *
   * Note that this formulation does not use the minimal 11 parameter 
   * representation. It also ignores the constraints on the rotation matrix
   * (R_3^T * R_3) = I, det(R_3) = 1. These are enforced after the computation.
   *
   * In addition, the use of the rotation matrix entries as our parameters 
   * requires us to extract the Euler angles. This is not a one to one mapping.
   * 
   * IN MOST CASES THERE ARE TWO POSSIBLE SETS OF EULER ANGLES THAT YIELD THE
   * SAME MATRIX DEPENDING ON THE SIGN OF THE SQRT IN THE COMPUTATION OF omega_y: 
   * 
   * omega_y = atan2(-R(2,0),+sqrt(R(0,0)^2 + R(1,0)^2)) or 
   *           atan2(-R(2,0), -sqrt(R(0,0)^2 + R(1,0)^2))
   * cy1 = cos(omega_y1)
   * cy2 = cos(omega_y2)
   * omega_x = atan2(R(2,1)/cy1, R(2,2)/cy1) or
   *           atan2(R(2,1)/cy2, R(2,2)/cy2)
   * omega_z = atan2(R(1,0)/cy1, R(0,0)/cy1) or
   *           atan2(R(1,0)/cy2, R(0,0)/cy2)
   *
   * when omega_y~pi/2 or omega_y~-pi/2 we arbitrarily set 
   * omega_z=0
   * omega_x = atan2(R3(0,1), R3(1,1))
   *
	 * @param data The US calibration should minimize the least squares 
   *             error of this data set.
	 * @param parameters This vector is cleared and then filled with the computed 
   *                   US calibration transformation parameters and the 
   *                   estimated 3D point location in the
   *                   world/tracker coordinate system (t_1).
   * [t_1x, t_1y, t_1z, t_3x, t_3y, t_3z, omega_z, omega_y, omega_x, m_x, m_y, m_x*R3(:,1), m_y*R3(:,2), R3(:,3)].
   *
   *                   Angles are in radians.
   *
   *                   The reason for providing the composition of rotation and 
   *                   scale, last nine entries, is computational efficiency of
   *                   the agree method. 
   *
   *                   If the vector contains less than four data elements or 
   *                   they do not provide a set of linearly independent equation
   *                   systems the resulting parameters vector is empty 
   *                   (size == 0). A natural example where this happens is when
   *                   the US probe's DRF is aligned with the tracking system's
   *                   coordinate system (R2i=I) and the probe is only 
   *                   translated, that is A=[ui*I, vi*I, I, -I], clearly 
   *                   rank(A)<=9.
	 */
	void analyticLeastSquaresEstimate(std::vector< DataType *> &data, 
                                    std::vector<double> &parameters);

  /**
	 * Compute a least squares estimate of the US calibration and target point
   * location in the world/tracker coordinate system, as defined by the given 
   * data elements (pairs of US-DRF transformation and 2D point location).
	 * This implementation is of an iterative least squares solution:
	 *    min (sum ( (A_11i*m_x*cz*cy +
   *                 A_12i*m_x*sz*cy -
   *                 A_13i*m_x*sy +
   *                 A_14i*m_y*(cz*sy*sx - sz*cx) +
   *                 A_15i*m_y*(sz*sy*sx+cz*cx) +
   *                 A_16i*m_y*cy*sx +
   *                 A_17i*t_3x +
   *                 A_18i*t_3y +
   *                 A_19i*t_3z -
   *                 t_1x +
   *                 t_2ix)^2 +
   *                 (A_21i*m_x*cz*cy +
   *                  A_22i*m_x*sz*cy -
   *                  A_23i*m_x*sy +
   *                  A_24i*m_y*(cz*sy*sx - sz*cx) +
   *                  A_25i*m_y*(sz*sy*sx+cz*cx) +
   *                  A_26i*m_y*cy*sx +
   *                  A_27i*t_3x +
   *                  A_28i*t_3y +
   *                  A_29i*t_3z -
   *                  t_1y +
   *                  t_2iy)^2 +
   *                  (A_31i*m_x*cz*cy +
   *                   A_32i*m_x*sz*cy -
   *                   A_33i*m_x*sy +
   *                   A_34i*m_y*(cz*sy*sx - sz*cx) +
   *                   A_35i*m_y*(sz*sy*sx+cz*cx) +
   *                   A_36i*m_y*cy*sx +
   *                   A_37i*t_3x +
   *                   A_38i*t_3y +
   *                   A_39i*t_3z -
   *                   t_1z +
   *                   t_2iz)^2 ) )
	 *
   * where cx - cos(omega_x), sx - sin(omega_x), A = [u_i*R_2i, v_i*R_2i, R_2i -I].
   * We use a minimal parameterization (11 parameters): 
   * [t_1 (x,y,z), t_3 (x,y,z), omega (z,y,x), m (x,y)]
   *
   * The Levenberg-Marquardt algorithm is used for the minimization.
   *
	 * @param data The US calibration should minimize the least squares 
   *             error of this data set.
	 * @param initialParameters This vector contains the initial parameter values 
   *                          for nonlinear estimation.
	 * @param finalParameters This vector is cleared and then filled with the computed 
   *                   US calibration transformation parameters and the 
   *                   estimated 3D point location in the
   *                   world/tracker coordinate system (t_1).
   * [t_1x, t_1y, t_1z, t_3x, t_3y, t_3z, omega_z, omega_y, omega_x, m_x, m_y, m_x*R3(:,1), m_y*R3(:,2), R3(:,3)].
   *
   *                   Angles are in radians.
   *
   *                   The reason for providing the composition of rotation and 
   *                   scale, last nine entries, is computational efficiency of
   *                   the agree method. 
   *
   * NOTE: As this method encapsulates nonlinear optimization (Levenberg-Marquardt) 
   *       it may not converge due to the specific choice of LM parameter settings.
   *       Possibly modify these values to be appropriate for user's specific input.
	 */
  void iterativeLeastSquaresEstimate(std::vector< DataType *> &data, 
																		 std::vector<double> &initialParameters, 
                                     std::vector<double> &finalParameters);

  static void getDistanceStatistics(const std::vector<double> &parameters, 
                                    const std::vector<DataType> &data, 
                                    std::vector<double> &distances,
                                    double &min, double &max, double &mean);

private:

  class SumSquaresCalibrationPointsDistanceFunction : 
   public vnl_least_squares_function {
		public:
     SumSquaresCalibrationPointsDistanceFunction(std::vector<DataType *> *data);
			void setData(std::vector<DataType *> *data);
      virtual void f( vnl_vector<double> const &x, vnl_vector<double> &fx );
      virtual void gradf( vnl_vector<double> const& x, vnl_matrix<double>& jacobian );
		private:
			std::vector<DataType *> *data;
	};

	double deltaSquared; 
           //analytic or iterative least squares
  LeastSquaresType lsType; 
};













/**
 * The data elements used by the pointer based estimator are triplets that 
 * include:
 * (a). T2 - the transformation mapping points from the US dynamic reference 
 *           frame to the fixed, tracker, coordinate system T.
 * (b). q - the target point's coordinates in the image coordinate system 
 *          [u,v,0].
 * (c). p - the target point's coordinates in the fixed, tracker, coordinate
 *          system
 */
struct CalibratedPointerTargetUSCalibrationParametersEstimatorDataType {
    Frame T2;  
    Point2D q;
    Point3D p;
};


class CalibratedPointerTargetUSCalibrationParametersEstimator : 
  public ParametersEstimator<CalibratedPointerTargetUSCalibrationParametersEstimatorDataType,double> {
public:
  
  enum LeastSquaresType {ANALYTIC = 0, ITERATIVE};

  typedef CalibratedPointerTargetUSCalibrationParametersEstimatorDataType DataType;

  /**
	 * Object constructor.
	 * @param delta A data element (triplet of: US-DRF transformation, 2D point 
   *              coordinates and 3D point coordinates) agrees with specific 
   *              calibration parameter values if the 3D distance between the 
   *              mapped point and its known location is less than 'delta'.
	 * @param lsType When the leastSquaresEstimate() method is called it computes 
   *               an analytic or iterative fit. This flag tells it which one.
   */
	CalibratedPointerTargetUSCalibrationParametersEstimator(double delta, 
                                                          LeastSquaresType lsType = ITERATIVE);

	/**
	 * Compute the affine US calibration transformation defined by the given data 
   * elements. 
   *
   * @param data A vector containing three data elements.
	 * @param parameters This vector is cleared and then filled with the computed 
   *                   transformation parameters 
   * [t_3x, t_3y, t_3z, omega_z, omega_y, omega_x, m_x, m_y, m_x*R3(:,1), m_y*R3(:,2), R3(:,3)].
   *
   *                   Angles are in radians.
   *
   *                   The reason for providing the composition of rotation and 
   *                   scale, last nine entries, is computational efficiency of
   *                   the agree method. 
   *
   *                   If the vector contains less than three data elements or 
   *                   they do not provide a set of linearly independent equation
   *                   systems the resulting parameters vector is empty 
   *                   (size == 0).
	 */	
  virtual void estimate(std::vector<DataType *> &data, 
                        std::vector<double> &parameters);
	virtual void estimate(std::vector<DataType> &data, 
                        std::vector<double> &parameters);


	/**
	 * Compute the least squares affine US calibration transformation defined by 
   * the given data elements. This may be either an analytic or iterative least 
   * squares computation depending on the LeastSquaresType settings. For further 
   * information see iterativeLeastSquaresEstimate(), and 
   * analyticLeastSquaresEstimate().
   *
   * @param data A vector containing at least three data elements.
	 * @param parameters This vector is cleared and then filled with the computed 
   *                   transformation parameters 
   * [t_3x, t_3y, t_3z, omega_z, omega_y, omega_x, m_x, m_y, m_x*R3(:,1), m_y*R3(:,2), R3(:,3)].
   *
   *                   Angles are in radians.
   *
   *                   The reason for providing the composition of rotation and 
   *                   scale, last nine entries, is computational efficiency of
   *                   the agree method. 
   *
   *                   If the vector contains less than four data elements or 
   *                   they do not provide a set of linearly independent equation
   *                   systems the resulting parameters vector is empty 
   *                   (size == 0).
	 */	
	virtual void leastSquaresEstimate(std::vector<DataType *> &data, 
                                    std::vector<double> &parameters);
	virtual void leastSquaresEstimate(std::vector<DataType> &data, 
                                    std::vector<double> &parameters);

  /**
   * A given data element is consistent, agrees, with the parameter values if 
   * the distance between the mapped 3D point and its known location is less 
   * than a user set threshold (see setDelta()).
   *
	 * @param parameters This vector contains the US calibration transformation 
   *                   parameters.
   * [t_3x, t_3y, t_3z, omega_z, omega_y, omega_x, m_x, m_y, m_x*R3(:,1), m_y*R3(:,2), R3(:,3)].
   *
   *                   Angles are in radians.
   *
   * @param data A data element comprised of the US-DRF transformation, the
   *             2D point location in the image coordinate system, and the 3D
   *             pointer tip location in the tracker coordinate system.
   *
   * Note: If the parameters vector is too short the method will throw an 
   *       exception as it tries to access a non existing entry. 
   *
   */
	virtual bool agree(std::vector<double> &parameters, DataType &data);

  /**
	 * Change the type of least squares solution.
	 * @param lsType When the leastSquaresEstimate() method is called it computes 
   *               an analytic or iterative fit. This flag tells it which one.
   */
	void setLeastSquaresType(LeastSquaresType lsType) {this->lsType = lsType;}

  /**
   * Set the threshold that defines if a data element agrees with a specific 
   * model estimate.
	 * @param delta A data element (triplet of US-DRF transformation, 2D image
   *              point location, and 3D pointer tip location) 
   *              agrees with specific calibration parameter values if
   *              the 3D distance between the mapped point and its known location
   *              is less than 'delta'.
   */
	void setDelta(double delta) {this->deltaSquared = delta*delta;}

	/**
	 * Compute a least squares estimate of the US calibration, as defined by the 
   * given data elements (triplets of US-DRF transformation, 2D point location,
   * and 3D pointer tip location).
	 * This implementation is of an analytic least squares error:
	 * min||Ax-b||, where A = [u_i*R_2i, v_i*R_2i, R_2i], 
   *                    x = [m_x*R_3(:,1), m_y*R_3(:,2), t_3], 
   *                    b = [p_i-t_2i]
	 *
   * Note that this formulation does not use the minimal 8 parameter 
   * representation. It also ignores the constraints on the rotation matrix
   * (R_3^T * R_3) = I, det(R_3) = 1. These are enforced after the computation.
   *
   * In addition, the use of the rotation matrix entries as our parameters 
   * requires us to extract the Euler angles. This is not a one to one mapping.
   * 
   * IN MOST CASES THERE ARE TWO POSSIBLE SETS OF EULER ANGLES THAT YIELD THE
   * SAME MATRIX DEPENDING ON THE SIGN OF THE SQRT IN THE COMPUTATION OF omega_y: 
   * 
   * omega_y = atan2(-R(2,0),+sqrt(R(0,0)^2 + R(1,0)^2)) or 
   *           atan2(-R(2,0), -sqrt(R(0,0)^2 + R(1,0)^2))
   * cy1 = cos(omega_y1)
   * cy2 = cos(omega_y2)
   * omega_x = atan2(R(2,1)/cy1, R(2,2)/cy1) or
   *           atan2(R(2,1)/cy2, R(2,2)/cy2)
   * omega_z = atan2(R(1,0)/cy1, R(0,0)/cy1) or
   *           atan2(R(1,0)/cy2, R(0,0)/cy2)
   *
   * when omega_y~pi/2 or omega_y~-pi/2 we arbitrarily set 
   * omega_z=0
   * omega_x = atan2(R3(0,1), R3(1,1))
   *
	 * @param data The US calibration should minimize the least squares 
   *             error of this data set.
	 * @param parameters This vector is cleared and then filled with the computed 
   *                   US calibration transformation parameters.
   * [t_3x, t_3y, t_3z, omega_z, omega_y, omega_x, m_x, m_y, m_x*R3(:,1), m_y*R3(:,2), R3(:,3)].
   *
   *                   Angles are in radians.
   *
   *                   The reason for providing the composition of rotation and 
   *                   scale, last nine entries, is computational efficiency of
   *                   the agree() method. 
   *
   *                   If the vector contains less than three data elements or 
   *                   they do not provide a set of linearly independent equation
   *                   systems the resulting parameters vector is empty 
   *                   (size == 0).
	 */
	void analyticLeastSquaresEstimate(std::vector< DataType *> &data, 
                                    std::vector<double> &parameters);

  /**
	 * Compute a least squares estimate of the US calibration and target point
   * location in the world/tracker coordinate system, as defined by the given 
   * data elements (triplets of US-DRF transformation, 2D point location and 
   * 3D pointer tip location).
	 * This implementation is of an iterative least squares solution:
	 *    min (sum ( (A_11i*m_x*cz*cy +
   *                 A_12i*m_x*sz*cy -
   *                 A_13i*m_x*sy +
   *                 A_14i*m_y*(cz*sy*sx - sz*cx) +
   *                 A_15i*m_y*(sz*sy*sx+cz*cx) +
   *                 A_16i*m_y*cy*sx +
   *                 A_17i*t_3x +
   *                 A_18i*t_3y +
   *                 A_19i*t_3z -
   *                 p_ix +
   *                 t_2ix)^2 +
   *                 (A_21i*m_x*cz*cy +
   *                  A_22i*m_x*sz*cy -
   *                  A_23i*m_x*sy +
   *                  A_24i*m_y*(cz*sy*sx - sz*cx) +
   *                  A_25i*m_y*(sz*sy*sx+cz*cx) +
   *                  A_26i*m_y*cy*sx +
   *                  A_27i*t_3x +
   *                  A_28i*t_3y +
   *                  A_29i*t_3z -
   *                  p_iy +
   *                  t_2iy)^2 +
   *                  (A_31i*m_x*cz*cy +
   *                   A_32i*m_x*sz*cy -
   *                   A_33i*m_x*sy +
   *                   A_34i*m_y*(cz*sy*sx - sz*cx) +
   *                   A_35i*m_y*(sz*sy*sx+cz*cx) +
   *                   A_36i*m_y*cy*sx +
   *                   A_37i*t_3x +
   *                   A_38i*t_3y +
   *                   A_39i*t_3z -
   *                   p_iz +
   *                   t_2iz)^2 ) )
	 *
   * where cx - cos(omega_x), sx - sin(omega_x), A = [u_i*R_2i, v_i*R_2i, R_2i].
   * We use a minimal parameterization (8 parameters): 
   * [t_3 (x,y,z), omega (z,y,x), m (x,y)].
   *
   * The Levenberg-Marquardt algorithm is used for the minimization.
   *
	 * @param data The US calibration should minimize the least squares 
   *             error of this data set.
	 * @param initialParameters This vector contains the initial parameter values 
   *                          for nonlinear estimation.
	 * @param finalParameters This vector is cleared and then filled with the computed 
   *                   US calibration transformation parameters.
   * [t_3x, t_3y, t_3z, omega_z, omega_y, omega_x, m_x, m_y, m_x*R3(:,1), m_y*R3(:,2), R3(:,3)].
   *
   *                   Angles are in radians.
   *
   *                   The reason for providing the composition of rotation and 
   *                   scale, last nine entries, is computational efficiency of
   *                   the agree method. 
   *
   * NOTE: As this method encapsulates nonlinear optimization (Levenberg-Marquardt) 
   *       it may not converge due to the specific choice of LM parameter settings.
   *       Possibly modify these values to be appropriate for user's specific input.
	 */
  void iterativeLeastSquaresEstimate(std::vector< DataType *> &data, 
																		 std::vector<double> &initialParameters, 
                                     std::vector<double> &finalParameters);

  static void getDistanceStatistics(const std::vector<double> &parameters, 
                                    const std::vector<DataType> &data, 
                                    std::vector<double> &distances,
                                    double &min, double &max, double &mean);

private:

  class SumSquaresCalibrationPointsDistanceFunction : 
   public vnl_least_squares_function {
		public:
     SumSquaresCalibrationPointsDistanceFunction(std::vector<DataType *> *data);
			void setData(std::vector<DataType *> *data);
      virtual void f( vnl_vector<double> const &x, vnl_vector<double> &fx );
      virtual void gradf( vnl_vector<double> const& x, vnl_matrix<double>& jacobian );
		private:
			std::vector<DataType *> *data;
	};

	double deltaSquared; 
           //analytic or iterative least squares
  LeastSquaresType lsType; 
};


} //namespace lsqrRecipes

#endif //_SINGLE_POINT_TARGET_US_CALIBRATION_PARAMETERS_ESTIMATOR_H_
