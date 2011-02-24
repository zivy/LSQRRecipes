#ifndef _PLANE_PHANTOM_US_CALIBRATION_PARAMETERS_ESTIMATOR_H_
#define _PLANE_PHANTOM_US_CALIBRATION_PARAMETERS_ESTIMATOR_H_

#include "copyright.h"

#include <vnl/vnl_least_squares_function.h>
#include "ParametersEstimator.h"
#include "Point2D.h"
#include "Frame.h"


/**
 * This class estimates the transformation mapping points from the US image
 * coordinate system to the coordinate system of a tracked reference frame that
 * is fixed relative to the US probe. The method is based on the use of a planar
 * phantom, either the bottom of a water bath or a planar membrane.
 *For additional information with regard to US calibration see: 
 *
 * "A review of calibration techniques for freehand 3-D ultrasound systems", 
 * L. Mercier, T. Lango, F. Lindseth, L.D. Collins, Ultrasound in Med. & Biol., 
 * vol. 31(2), pp. 143-165, 2005.
 * 
 * @author: Ziv Yaniv (zivy@isis.georgetown.edu)
 *
 */

namespace lsqrRecipes {

/**
 * The data elements used by the plane based estimator are pairs that 
 * include:
 * (a). T2 - the transformation mapping points from the US dynamic reference 
 *           frame to the fixed, tracker, coordinate system T.
 * (b). q - 2D coordinates in the image coordinate system [u,v,0] corresponding
 *          to a 3D point on the imaged plane.
 */
struct PlanePhantomUSCalibrationParametersEstimatorDataType {
    Frame T2;  
    Point2D q;
};


class PlanePhantomUSCalibrationParametersEstimator : 
  public ParametersEstimator<PlanePhantomUSCalibrationParametersEstimatorDataType,double> {
public:
  
  enum LeastSquaresType {ANALYTIC = 0, ITERATIVE};

  typedef PlanePhantomUSCalibrationParametersEstimatorDataType DataType;

  /**
	 * Object constructor.
	 * @param delta A data element (pair of US-DRF transformation and 2D point location) 
   *              agrees with specific calibration parameter values if
   *              the distance between the mapped point and the x-y plane (z=0)
   *              is less than 'delta'.
	 * @param lsType When the leastSquaresEstimate() method is called it computes 
   *               an analytic or iterative fit. This flag tells it which one.
   */
	PlanePhantomUSCalibrationParametersEstimator(double delta, 
                                               LeastSquaresType lsType = ITERATIVE);

	/**
	 * Compute the affine US calibration transformation defined by the given data 
   * elements. 
   *
   * @param data A vector containing 31 data elements.
	 * @param parameters This vector is cleared and then filled with the computed 
   *                   transformation parameters 
   * [omega1_y, omega1_x, t_1z, t_3x, t_3y, t_3z, omega3_z, omega3_y, omega3_x, m_x, m_y, [*]].
   * [*] - 30 entries dependent on the previous 11:
   * m_x*R1(3,1)*R3(:,1),
   * m_x*R1(3,2)*R3(:,1),
   * m_x*R1(3,3)*R3(:,1),
   * m_y*R1(3,1)*R3(:,2),
   * m_y*R1(3,2)*R3(:,2),
   * m_y*R1(3,3)*R3(:,2),
   * R1(3,1)*t_3,
   * R1(3,2)*t_3,
   * R1(3,3)*t_3,
   * R1(3,:)
   * 
   *                   Angles are in radians.
   *
   *                   The reason for providing the last 30 entries, is 
   *                   computational efficiency of the agree method. 
   *
   *                   If the vector contains less than 31 data elements or 
   *                   they do not provide a set of linearly independent equation
   *                   systems the resulting parameters vector is empty 
   *                   (size == 0).
	 */	
  virtual void estimate(std::vector<DataType *> &data, 
                        std::vector<double> &parameters);
	virtual void estimate(std::vector<DataType> &data, 
                        std::vector<double> &parameters);


	virtual void leastSquaresEstimate(std::vector<DataType *> &data, 
                                    std::vector<double> &parameters);
	virtual void leastSquaresEstimate(std::vector<DataType> &data, 
                                    std::vector<double> &parameters);

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
   *              the distance between the mapped point and the x-y plane (z=0)
   *              is less than 'delta'.
   */
	void setDelta(double delta) {this->deltaSquared = delta*delta;}

	void analyticLeastSquaresEstimate(std::vector< DataType *> &data, 
                                    std::vector<double> &parameters);

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

#endif //_PLANE_PHANTOM_US_CALIBRATION_PARAMETERS_ESTIMATOR_H_
