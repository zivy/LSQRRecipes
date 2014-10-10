#ifndef _PIVOT_CALIBRATION_PARAMETERS_ESTIMATOR_H_
#define _PIVOT_CALIBRATION_PARAMETERS_ESTIMATOR_H_

#include "copyright.h"

#include "ParametersEstimator.h"
#include "Frame.h"

/**
 * This class performs pivot calibration using an algebraic method which only
 * requires solving a single set of linear equations. It is based on the 
 * observation that when we are pivoting around a fixed point all transformations 
 * yield the following three equations:
 * R_i*{DRF}^t + t_i = W^t
 *
 * where [R_i,t_i] is the transformation from the tracking system to the Dynamic
 * Reference Frame (DRF), {DRF}^t is the translation from the DRF to the pivoting
 * point and W^t is the translation from the tracking system to the pivoting point.
 * To obtain {DRF}^t and W^t we solve the overdetermined equation system:
 * [R_1 -I]  [{DRF}^t]   [-t_1]
 * [   :  ]  [       ] = [ :  ]
 * [R_m -I]  [  W^t  ]   [-t_m]
 *
 * For additional details see:
 * Z. Yaniv, "Which Pivot Calibration?", SPIE Medical Imaging: Image-Guided 
 * Procedures, Robotic Interventions, and Modeling", 2015.
 *
 * @author: Ziv Yaniv 
 */
namespace lsqrRecipes {

class PivotCalibrationEstimator : public ParametersEstimator< Frame ,double> {
public:

  /**
	 * Object constructor.
	 * @param delta A transformation is consistent with a given estimate of 
   *              [{DRF}^t, W^t] if ||R*{DRF}^t + t - W^t||<delta.
   */
  PivotCalibrationEstimator(double delta);

	/**
	 * Compute the two translation vectors defined by the given transformations.
	 * @param data A vector containing three transformations.
	 * @param parameters This vector is cleared and then filled with the values,
   *                   [{DRF}^t_x, {DRF}^t_y, {DRF}^t_z, W^t_x, W^t_y, W^t_z]). 
   *                   If the data contains less than three 
   *                   transformations or the equation system that they define is
   *                   singular then the parameters vector is empty (size = 0).
	 */
	virtual void estimate(std::vector< Frame *> &data, 
                        std::vector<double> &parameters);
	virtual void estimate(std::vector< Frame > &data, 
                        std::vector<double> &parameters);

	/**
	 * Compute a least squares estimate of the two translation vectors defined by 
   * the given transformations.
	 * @param data A vector containing three or more transformations.
	 * @param parameters This vector is cleared and then filled with the values,
   *                   [{DRF}^t_x, {DRF}^t_y, {DRF}^t_z, W^t_x, W^t_y, W^t_z]).
   *                   If the data contains less than three transformations or 
   *                   the equation system that they define is singular then the 
   *                   parameters vector is empty (size = 0).
	 */
	virtual void leastSquaresEstimate(std::vector< Frame *> &data, 
                                    std::vector<double> &parameters);
	virtual void leastSquaresEstimate(std::vector< Frame > &data, 
                                    std::vector<double> &parameters);

	/**
   * Return true if the given transformation is consistent with the translations 
   * defined by the parameters, error is smaller than 'delta' (see constructor).
   *
	 * @param parameters The translations  [{DRF}^t, W^t].
	 * @param data Check that the distance between this point and the hyperplane is 
   *             less than 'delta'.
	 * @param delta A transformation is consistent with a given estimate of 
   *              [{DRF}^t, W^t] if ||R*{DRF}^t + t - W^t||<delta.
   */
	virtual bool agree(std::vector<double> &parameters, Frame &data);

	/**
	 * Change the error defining if a transformation is consistent with the given
   * parameter values.
	 * @param delta A transformation is consistent with a given estimate of 
   *              [{DRF}^t, W^t] if ||R*{DRF}^t + t - W^t||<delta.
	 */	void setDelta(double delta) {this->delta = delta;}

private:
    //given a transformation/frame and the two translations estimated by pivot
    //calibration the estimated translations are consistent with the 
    //transformation if || R{DRF}^t + t - W^t ||<delta
	double delta; 
};

} //namespace lsqrRecipes


#endif //_PIVOT_CALIBRATION_PARAMETERS_ESTIMATOR_H_
