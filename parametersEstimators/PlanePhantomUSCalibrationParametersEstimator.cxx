#include "PlanePhantomUSCalibrationParametersEstimator.h"
#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_cross.h>
#include <vnl/algo/vnl_levenberg_marquardt.h>
#include "Epsilon.h"

namespace lsqrRecipes {

PlanePhantomUSCalibrationParametersEstimator::PlanePhantomUSCalibrationParametersEstimator(double delta, LeastSquaresType lsType) : 
  ParametersEstimator<DataType,double>(31) 
{
  this->deltaSquared = delta*delta;
  this->lsType = lsType;
}
/*****************************************************************************/
void PlanePhantomUSCalibrationParametersEstimator::estimate(std::vector<DataType *> &data, 
																	                          std::vector<double> &parameters)
{
  parameters.clear();
  if(data.size() != this->minForEstimate)
    return;
        //forward the work to the analytic least squares
  analyticLeastSquaresEstimate(data,parameters);
}
/*****************************************************************************/
void PlanePhantomUSCalibrationParametersEstimator::estimate(std::vector<DataType> &data, 
																	                          std::vector<double> &parameters)
{
	std::vector<DataType *> usedData;
	size_t dataSize = data.size();
	for(size_t i=0; i<dataSize; i++)
		usedData.push_back(&(data[i]));
	estimate(usedData,parameters);
}
/*****************************************************************************/
void PlanePhantomUSCalibrationParametersEstimator::leastSquaresEstimate(std::vector<DataType *> &data, 
																							                          std::vector<double> &parameters)
{
	parameters.clear();
	if(data.size()<this->minForEstimate)
		return;

  std::vector<double> initialParameters;
	switch(this->lsType) {
		case ANALYTIC: 
			analyticLeastSquaresEstimate(data, parameters);
			break;
		case ITERATIVE:
	                   //analytic least squares for initial estimate         
	    analyticLeastSquaresEstimate(data, initialParameters);
	    if(initialParameters.size() == 0)
		    return;
			iterativeLeastSquaresEstimate(data, initialParameters, parameters);
			break;
	}

}
/*****************************************************************************/
void PlanePhantomUSCalibrationParametersEstimator::leastSquaresEstimate(std::vector<DataType> &data, 
																							                          std::vector<double> &parameters)
{
	std::vector<DataType *> usedData;
	size_t dataSize = data.size();
	for(size_t i=0; i<dataSize; i++)
		usedData.push_back(&(data[i]));
	leastSquaresEstimate(usedData,parameters);
}
/*****************************************************************************/
/*
 * If the parameters vector is too short an exception will
 * be thrown by the vector's [] operator. We don't check vector length explicitly.
 */
bool PlanePhantomUSCalibrationParametersEstimator::agree(std::vector<double> &parameters, 
                                                         DataType &data)
{
  double R2[3][3], t2[3], u, v;
  double err;

  data.T2.getRotationMatrix(R2);
  data.T2.getTranslation(t2);
  u = data.q[0];
  v = data.q[1];

  err = u*R2[0][0]*parameters[11] + 
        u*R2[0][1]*parameters[12] + 
        u*R2[0][2]*parameters[13] +
        u*R2[1][0]*parameters[14] + 
        u*R2[1][1]*parameters[15] + 
        u*R2[1][2]*parameters[16] +
        u*R2[2][0]*parameters[17] + 
        u*R2[2][1]*parameters[18] + 
        u*R2[2][2]*parameters[19] +
        v*R2[0][0]*parameters[20] + 
        v*R2[0][1]*parameters[21] + 
        v*R2[0][2]*parameters[22] +
        v*R2[1][0]*parameters[23] + 
        v*R2[1][1]*parameters[24] + 
        v*R2[1][2]*parameters[25] +
        v*R2[2][0]*parameters[26] + 
        v*R2[2][1]*parameters[27] + 
        v*R2[2][2]*parameters[28] +
        R2[0][0]*parameters[29] + 
        R2[0][1]*parameters[30] + 
        R2[0][2]*parameters[31] +
        R2[1][0]*parameters[32] + 
        R2[1][1]*parameters[33] + 
        R2[1][2]*parameters[34] +
        R2[2][0]*parameters[35] + 
        R2[2][1]*parameters[36] + 
        R2[2][2]*parameters[37] +
        t2[0]*parameters[38] +
        t2[1]*parameters[39] +
        t2[2]*parameters[40] +
        parameters[2];

  return (err*err<this->deltaSquared);
}
/*****************************************************************************/
/**
 * Solve the homogenous equation system Ax=0, where A is a 3nX31 matrix (n>=31), 
 * and x is a 31x1. Each data element provides one equation as follows:
 * 
 *
 * [ u_i*R2_i(1,:) ] ^T  [m_x*R1(3,1)*R3(:,1)]  
 * [ u_i*R2_i(2,:) ]     [m_x*R1(3,2)*R3(:,1)]  
 * [ u_i*R2_i(3,:) ]     [m_x*R1(3,3)*R3(:,1)]   
 * [ v_i*R2_i(1,:) ]     [m_y*R1(3,1)*R3(:,2)]  
 * [ v_i*R2_i(2,:) ]     [m_y*R1(3,2)*R3(:,2)]  
 * [ v_i*R2_i(3,:) ]     [m_y*R1(3,3)*R3(:,2)]   = [0]  
 * [   R2_i(1,:)   ]     [   R1(3,1)*t3(:)   ]  
 * [   R2_i(2,:)   ]     [   R1(3,2)*t3(:)   ]  
 * [   R2_i(3,:)   ]     [   R1(3,3)*t3(:)   ]  
 * [   t2(:)       ]     [   R1(3,:)         ]
 * [     1         ]     [     t1_z          ]
 *
 */
void PlanePhantomUSCalibrationParametersEstimator::analyticLeastSquaresEstimate(std::vector< DataType *> &data, 
                                                                                std::vector<double> &parameters)
{
	parameters.clear();
           //not enough data elements for estimate
	if(data.size()<this->minForEstimate) 
		return;

  unsigned int i, numPoints;
  double R1_31, R1_32, R1_33, omega1_x, omega1_y, t1_z;
  double omega3_x, omega3_y, omega3_z, m_x, m_y, t3_x, t3_y, t3_z;

  numPoints = static_cast<unsigned int>(data.size());
  vnl_matrix<double> A(numPoints,31);
  vnl_vector<double> x(31);
  double R2[3][3], t2x, t2y, t2z, ui, vi;

  for(i=0; i<numPoints; i++) {
           //transformation from US reference frame to tracker
    lsqrRecipes::Frame &frm = data[i]->T2;               
    frm.getRotationMatrix(R2);
    frm.getTranslation(t2x,t2y,t2z);
           //2D pixel coordinates
    ui = data[i]->q[0];
    vi = data[i]->q[1];
           //create row in the matrix A
    A(i,0) = R2[0][0]*ui; 
    A(i,1) = R2[0][1]*ui;
    A(i,2) = R2[0][2]*ui;
    A(i,3) = R2[1][0]*ui; 
    A(i,4) = R2[1][1]*ui; 
    A(i,5) = R2[1][2]*ui; 
    A(i,6) = R2[2][0]*ui; 
    A(i,7) = R2[2][1]*ui; 
    A(i,8) = R2[2][2]*ui;
    A(i,9) = R2[0][0]*vi; 
    A(i,10) = R2[0][1]*vi;
    A(i,11) = R2[0][2]*vi;
    A(i,12) = R2[1][0]*vi; 
    A(i,13) = R2[1][1]*vi; 
    A(i,14) = R2[1][2]*vi; 
    A(i,15) = R2[2][0]*vi; 
    A(i,16) = R2[2][1]*vi; 
    A(i,17) = R2[2][2]*vi;
    A(i,18) = R2[0][0]; 
    A(i,19) = R2[0][1];
    A(i,20) = R2[0][2];
    A(i,21) = R2[1][0]; 
    A(i,22) = R2[1][1]; 
    A(i,23) = R2[1][2]; 
    A(i,24) = R2[2][0]; 
    A(i,25) = R2[2][1]; 
    A(i,26) = R2[2][2];
    A(i,27) = t2x; 
    A(i,28) = t2y; 
    A(i,29) = t2z;
    A(i,30) = 1;    
  }

  vnl_svd<double> svdA(A);
  if(svdA.rank()<31) //points do not yield a solution
	  return;
           //get the right singular vector corresponding to the smallest singular
           //value, I assume the vectors in V correspond to a descending order
           //of the singular values (i.e. last column of V corresponds to smallest 
           //singular value).
  x = svdA.V().get_column(30);

  std::cout<<x<<"\n";

           //extract the calibration parameters:
           //we have the solution to Ax=0 s.t. ||x||=1. the solution we seek
           //also requires that ||R1(3,:)|| = ||x(27)*x(27)+x(28)*x(28)+x(29)*x(29)|| = 1
           //if we scale x accordingly we still have a valid solution with the
           //correct norm for the last row of the rotation matrix

           //the solution to the equation system can be invalid as we did not
           //enforce the constraints in its formulation
  double denominator = sqrt(x(27)*x(27)+x(28)*x(28)+x(29)*x(29));
  if(denominator<EPS)
    return;

           //NOTE: scaleFactor can be either +- 1/denominator, the values for
           //      T3 are invariant to this scaling, the values of T1 are not
           //      we arbitrarly choose to take the positive scale factor.
  double scaleFactor = 1/denominator;
  x*= scaleFactor;
  
  t1_z = x(30);  

  double smallAngle = 0.008726535498373935;//0.5 degrees
  double halfPI = 1.5707963267948966192313216916398;

          //last row of R1, get the two Euler rotation angles
  R1_31 = x(27);
  R1_32 = x(28);
  R1_33 = x(29);

           //two options for omega1_y depending on choice of +-1 for sqrt
           //we arbitrarily choose +1
  omega1_y = atan2(-R1_31, sqrt(R1_32*R1_32 + R1_33*R1_33));
             //omega1_x depends on the omega1_y value
	if(fabs(omega1_y - halfPI) > smallAngle && fabs(omega1_y + halfPI) > smallAngle) {
    double cy = cos(omega1_y);
    omega1_x = atan2(R1_32/cy, R1_33/cy);  
  }     //gimbal lock, omega1_y is approximatly plus/minus half PI, arbitrarily 
        //set omega1_x to 0 
  else {          
    omega1_x = 0.0;
  }

            //get the translation t3, "average out" the errors
  t3_x = (x(18)/R1_31 + x(21)/R1_32 + x(24)/R1_33)/3.0;
  t3_y = (x(19)/R1_31 + x(22)/R1_32 + x(25)/R1_33)/3.0;
  t3_z = (x(20)/R1_31 + x(23)/R1_32 + x(26)/R1_33)/3.0;

            //get the scale factors and the rotation matrix R3
  vnl_vector<double> r1a(3), r1b(3), r1c(3), r1(3), r2a(3), r2b(3), r2c(3), r2(3), r3(3);
  vnl_matrix<double> R3(3,3);

          //get the x scale factor and the first column of the rotation matrix 
  r1a(0) = x(0)/R1_31;
  r1a(1) = x(1)/R1_31;
  r1a(2) = x(2)/R1_31;
  r1b(0) = x(3)/R1_32;
  r1b(1) = x(4)/R1_32;
  r1b(2) = x(5)/R1_32;
  r1c(0) = x(6)/R1_33;
  r1c(1) = x(7)/R1_33;
  r1c(2) = x(8)/R1_33;
  r1 = (r1a+r1b+r1c)/3.0;
  m_x = r1.two_norm();
  r1.normalize();
          //get the y scale factor and the second column of the rotation matrix 
  r2a(0) = x(9)/R1_31;
  r2a(1) = x(10)/R1_31;
  r2a(2) = x(11)/R1_31;
  r2b(0) = x(12)/R1_32;
  r2b(1) = x(13)/R1_32;
  r2b(2) = x(14)/R1_32;
  r2c(0) = x(15)/R1_33;
  r2c(1) = x(16)/R1_33;
  r2c(2) = x(17)/R1_33;
  r2 = (r2a+r2b+r2c)/3.0;
  m_y = r2.two_norm();
  r2.normalize();
           //get the third column of the rotation matrix
  r3 = vnl_cross_3d(r1,r2);
          //the matrix R3=[r1,r2,r3] is not necessarily a rotation 
          //matrix, the orthonormality constraints were not enforced as part of 
          //the solution, get the closest (Frobenius norm) rotation matrix via
          //SVD
  R3.set_column(0,r1);
  R3.set_column(1,r2);
  R3.set_column(2,r3);
  vnl_svd<double> svdR3(R3);
  R3 = svdR3.U()*svdR3.V().transpose();
           //extract the Euler angles

            //two options for omega_y depending on choice of +-1 for sqrt
           //we arbitrarily choose +1
  omega3_y = atan2(-R3(2,0), sqrt(R3(0,0)*R3(0,0) + R3(1,0)*R3(1,0)));
             //omega_z and omega_x depend on the omega_y value
	if(fabs(omega3_y - halfPI) > smallAngle && fabs(omega3_y + halfPI) > smallAngle) {
    double cy = cos(omega3_y);
    omega3_z = atan2(R3(1,0)/cy, R3(0,0)/cy);
    omega3_x = atan2(R3(2,1)/cy, R3(2,2)/cy);
  }     //gimbal lock, omega_y is approximatly plus/minus half PI, arbitrarily 
        //set omega_z to 0 and compute omega_x
  else {          
    omega3_z = 0;
    omega3_x = atan2(R3(0,1), R3(1,1));
  }

       //set the parameters    
  parameters.push_back(omega1_y);
  parameters.push_back(omega1_x);
  parameters.push_back(t1_z); 
  parameters.push_back(t3_x);  
  parameters.push_back(t3_y);
  parameters.push_back(t3_z);
  parameters.push_back(omega3_z);  
  parameters.push_back(omega3_y);
  parameters.push_back(omega3_x);
  parameters.push_back(m_x);
  parameters.push_back(m_y);
                    //combinations of the minimal parameterization useful for
                    //repeated computations (i.e. the agree() method)
  parameters.push_back(m_x*R3(0,0)*R1_31);
  parameters.push_back(m_x*R3(1,0)*R1_31);
  parameters.push_back(m_x*R3(2,0)*R1_31);
  parameters.push_back(m_x*R3(0,0)*R1_32);
  parameters.push_back(m_x*R3(1,0)*R1_32);
  parameters.push_back(m_x*R3(2,0)*R1_32);
  parameters.push_back(m_x*R3(0,0)*R1_33);
  parameters.push_back(m_x*R3(1,0)*R1_33);
  parameters.push_back(m_x*R3(2,0)*R1_33);
  parameters.push_back(m_y*R3(0,1)*R1_31);
  parameters.push_back(m_y*R3(1,1)*R1_31);
  parameters.push_back(m_y*R3(2,1)*R1_31);
  parameters.push_back(m_y*R3(0,1)*R1_32);
  parameters.push_back(m_y*R3(1,1)*R1_32);
  parameters.push_back(m_y*R3(2,1)*R1_32);
  parameters.push_back(m_y*R3(0,1)*R1_33);
  parameters.push_back(m_y*R3(1,1)*R1_33);
  parameters.push_back(m_y*R3(2,1)*R1_33);
  parameters.push_back(t3_x*R1_31);  
  parameters.push_back(t3_y*R1_31);
  parameters.push_back(t3_z*R1_31);
  parameters.push_back(t3_x*R1_32);  
  parameters.push_back(t3_y*R1_32);
  parameters.push_back(t3_z*R1_32);
  parameters.push_back(t3_x*R1_33);  
  parameters.push_back(t3_y*R1_33);
  parameters.push_back(t3_z*R1_33);
  parameters.push_back(R1_31);  
  parameters.push_back(R1_32);
  parameters.push_back(R1_33);
}
/*****************************************************************************/
void PlanePhantomUSCalibrationParametersEstimator::iterativeLeastSquaresEstimate(std::vector< DataType *> &data, 
																	                                               std::vector<double> &initialParameters, 
                                                                                 std::vector<double> &finalParameters)
{
  unsigned int i, numParameters = 11;
	vnl_vector<double> parameters(numParameters);
	
	for(i=0; i<numParameters; i++)
	  parameters[i] = initialParameters[i];

	PlanePhantomUSCalibrationParametersEstimator::SumSquaresCalibrationPointsDistanceFunction 
    optimizedFunction(&data);		
  vnl_levenberg_marquardt lmOptimization(optimizedFunction);

            //set all the optimization tolerances (see vnl_nonlinear_minimizer)
	double gradTolerance = 10e-16;
	double parametersChangeTolerance = 10e-16;
  double functionChangeTolerance = 10e-16;
	int maxIterations = 5000;

  lmOptimization.set_f_tolerance(functionChangeTolerance);
  lmOptimization.set_x_tolerance(parametersChangeTolerance);
  lmOptimization.set_g_tolerance(gradTolerance);
  lmOptimization.set_max_function_evals(maxIterations);

	bool ok = lmOptimization.minimize(parameters);

  finalParameters.clear();
  if(ok) {
	  for(i=0; i<numParameters; i++)
	    finalParameters.push_back(parameters[i]);
                    //combinations of the minimal parameterization useful for
                    //repeated computations (i.e. the agree() method)
    double cx, cy, cz, sx, sy, sz;
    double R1_31, R1_32, R1_33, t3_x, t3_y, t3_z, m_x, m_y;  
    double R3_11, R3_12, R3_13, R3_21, R3_22, R3_23, R3_31, R3_32, R3_33;

    t3_x = finalParameters[3];
    t3_y = finalParameters[4];
    t3_z = finalParameters[5];
    m_x = finalParameters[9];
    m_y = finalParameters[10];

           //last row of R1
    cy = cos(finalParameters[0]);
    sy = sin(finalParameters[0]);
    cx = cos(finalParameters[1]);
    sx = sin(finalParameters[1]);
    R1_31 = -sy;
    R1_32 = cy*sx;
    R1_33 = cy*cx;

            //R3
    cz = cos(finalParameters[6]);
    sz = sin(finalParameters[6]);
    cy = cos(finalParameters[7]);
    sy = sin(finalParameters[7]);
    cx = cos(finalParameters[8]);
    sx = sin(finalParameters[8]);
    
	  R3_11 = cz*cy;    R3_12 = cz*sy*sx - sz*cx;    R3_13 = cz*sy*cx+sz*sx;
    R3_21 = sz*cy;    R3_22 = sz*sy*sx + cz*cx;    R3_23 = sz*sy*cx - cz*sx;
    R3_31 = -sy;      R3_32 = cy*sx;               R3_33 = cy*cx;

    finalParameters.push_back(m_x*R3_11*R1_31);
    finalParameters.push_back(m_x*R3_21*R1_31);
    finalParameters.push_back(m_x*R3_31*R1_31);
    finalParameters.push_back(m_x*R3_11*R1_32);
    finalParameters.push_back(m_x*R3_21*R1_32);
    finalParameters.push_back(m_x*R3_31*R1_32);
    finalParameters.push_back(m_x*R3_11*R1_33);
    finalParameters.push_back(m_x*R3_21*R1_33);
    finalParameters.push_back(m_x*R3_31*R1_33);
    finalParameters.push_back(m_y*R3_12*R1_31);
    finalParameters.push_back(m_y*R3_22*R1_31);
    finalParameters.push_back(m_y*R3_32*R1_31);
    finalParameters.push_back(m_y*R3_12*R1_32);
    finalParameters.push_back(m_y*R3_22*R1_32);
    finalParameters.push_back(m_y*R3_32*R1_32);
    finalParameters.push_back(m_y*R3_12*R1_33);
    finalParameters.push_back(m_y*R3_22*R1_33);
    finalParameters.push_back(m_y*R3_32*R1_33);
    finalParameters.push_back(t3_x*R1_31);  
    finalParameters.push_back(t3_y*R1_31);
    finalParameters.push_back(t3_z*R1_31);
    finalParameters.push_back(t3_x*R1_32);  
    finalParameters.push_back(t3_y*R1_32);
    finalParameters.push_back(t3_z*R1_32);
    finalParameters.push_back(t3_x*R1_33);  
    finalParameters.push_back(t3_y*R1_33);
    finalParameters.push_back(t3_z*R1_33);
    finalParameters.push_back(R1_31);  
    finalParameters.push_back(R1_32);
    finalParameters.push_back(R1_33);
  }

}
/*****************************************************************************/
void PlanePhantomUSCalibrationParametersEstimator::getDistanceStatistics(const std::vector<double> &parameters, 
                                                                         const std::vector<DataType> &data, 
                                                                         std::vector<double> &distances,
                                                                         double &min, double &max, double &mean)
{
  size_t n, i;
  double R2[3][3], t2[3], u, v;
  double err;

  data[0].T2.getRotationMatrix(R2);
  data[0].T2.getTranslation(t2);
  u = data[0].q[0];
  v = data[0].q[1];

  err = fabs(u*R2[0][0]*parameters[11] + 
        u*R2[0][1]*parameters[12] + 
        u*R2[0][2]*parameters[13] +
        u*R2[1][0]*parameters[14] + 
        u*R2[1][1]*parameters[15] + 
        u*R2[1][2]*parameters[16] +
        u*R2[2][0]*parameters[17] + 
        u*R2[2][1]*parameters[18] + 
        u*R2[2][2]*parameters[19] +
        v*R2[0][0]*parameters[20] + 
        v*R2[0][1]*parameters[21] + 
        v*R2[0][2]*parameters[22] +
        v*R2[1][0]*parameters[23] + 
        v*R2[1][1]*parameters[24] + 
        v*R2[1][2]*parameters[25] +
        v*R2[2][0]*parameters[26] + 
        v*R2[2][1]*parameters[27] + 
        v*R2[2][2]*parameters[28] +
        R2[0][0]*parameters[29] + 
        R2[0][1]*parameters[30] + 
        R2[0][2]*parameters[31] +
        R2[1][0]*parameters[32] + 
        R2[1][1]*parameters[33] + 
        R2[1][2]*parameters[34] +
        R2[2][0]*parameters[35] + 
        R2[2][1]*parameters[36] + 
        R2[2][2]*parameters[37] +
        t2[0]*parameters[38] +
        t2[1]*parameters[39] +
        t2[2]*parameters[40] +
        parameters[2]);

  distances.push_back(err);
  min = max = mean = distances[0];
  n = data.size();
  for(i=1; i<n; i++) {
    data[i].T2.getRotationMatrix(R2);
    data[i].T2.getTranslation(t2);
    u = data[i].q[0];
    v = data[i].q[1];

    err = fabs(u*R2[0][0]*parameters[11] + 
          u*R2[0][1]*parameters[12] + 
          u*R2[0][2]*parameters[13] +
          u*R2[1][0]*parameters[14] + 
          u*R2[1][1]*parameters[15] + 
          u*R2[1][2]*parameters[16] +
          u*R2[2][0]*parameters[17] + 
          u*R2[2][1]*parameters[18] + 
          u*R2[2][2]*parameters[19] +
          v*R2[0][0]*parameters[20] + 
          v*R2[0][1]*parameters[21] + 
          v*R2[0][2]*parameters[22] +
          v*R2[1][0]*parameters[23] + 
          v*R2[1][1]*parameters[24] + 
          v*R2[1][2]*parameters[25] +
          v*R2[2][0]*parameters[26] + 
          v*R2[2][1]*parameters[27] + 
          v*R2[2][2]*parameters[28] +
          R2[0][0]*parameters[29] + 
          R2[0][1]*parameters[30] + 
          R2[0][2]*parameters[31] +
          R2[1][0]*parameters[32] + 
          R2[1][1]*parameters[33] + 
          R2[1][2]*parameters[34] +
          R2[2][0]*parameters[35] + 
          R2[2][1]*parameters[36] + 
          R2[2][2]*parameters[37] +
          t2[0]*parameters[38] +
          t2[1]*parameters[39] +
          t2[2]*parameters[40] +
          parameters[2]);
    distances.push_back(err);    
    mean+=err;
    if(err>max)
      max = err;
    else if(err<min)
      min = err;
  }
  mean/= n;
}
/*****************************************************************************/
PlanePhantomUSCalibrationParametersEstimator::SumSquaresCalibrationPointsDistanceFunction::SumSquaresCalibrationPointsDistanceFunction(std::vector<DataType *> *data) :
vnl_least_squares_function(11,static_cast<unsigned int>(data->size()), vnl_least_squares_function::no_gradient/*vnl_least_squares_function::use_gradient*/)	                                                  
{
	this->data = data;
}

/*****************************************************************************/
void PlanePhantomUSCalibrationParametersEstimator::SumSquaresCalibrationPointsDistanceFunction::setData(std::vector<DataType *> *data)
{
	this->data = data;																											
}
/*****************************************************************************/
//[omega1_y, omega1_x, t1_x, t3_x, t3_y, t3_z, omega3_z, omega3_y, omage3_x, m_x, m_y]
void PlanePhantomUSCalibrationParametersEstimator::SumSquaresCalibrationPointsDistanceFunction::f(vnl_vector<double> const &x, vnl_vector<double> &fx)
{
    double cx, cy, cz, sx, sy, sz;
    double R1_31, R1_32, R1_33, t1_z, t3_x, t3_y, t3_z, m_x, m_y;  
    double R3_11, R3_12, R3_13, R3_21, R3_22, R3_23, R3_31, R3_32, R3_33;
    double R2[3][3], t2[3], u, v;

    t1_z = x[2];
    t3_x = x[3];
    t3_y = x[4];
    t3_z = x[5];
    m_x = x[9];
    m_y = x[10];

           //last row of R1
    cy = cos(x[0]);
    sy = sin(x[0]);
    cx = cos(x[1]);
    sx = sin(x[1]);
    R1_31 = -sy;
    R1_32 = cy*sx;
    R1_33 = cy*cx;

            //R3
    cz = cos(x[6]);
    sz = sin(x[6]);
    cy = cos(x[7]);
    sy = sin(x[7]);
    cx = cos(x[8]);
    sx = sin(x[8]);
    
	  R3_11 = cz*cy;    R3_12 = cz*sy*sx - sz*cx;    R3_13 = cz*sy*cx+sz*sx;
    R3_21 = sz*cy;    R3_22 = sz*sy*sx + cz*cx;    R3_23 = sz*sy*cx - cz*sx;
    R3_31 = -sy;      R3_32 = cy*sx;               R3_33 = cy*cx;

    double expr[27];

    expr[0] = m_x*R3_11*R1_31;
    expr[1] = m_x*R3_21*R1_31;
    expr[2] = m_x*R3_31*R1_31;
    expr[3] = m_x*R3_11*R1_32;
    expr[4] = m_x*R3_21*R1_32;
    expr[5] = m_x*R3_31*R1_32;
    expr[6] = m_x*R3_11*R1_33;
    expr[7] = m_x*R3_21*R1_33;
    expr[8] = m_x*R3_31*R1_33;
    expr[9] = m_y*R3_12*R1_31;
    expr[10] = m_y*R3_22*R1_31;
    expr[11] = m_y*R3_32*R1_31;
    expr[12] = m_y*R3_12*R1_32;
    expr[13] = m_y*R3_22*R1_32;
    expr[14] = m_y*R3_32*R1_32;
    expr[15] = m_y*R3_12*R1_33;
    expr[16] = m_y*R3_22*R1_33;
    expr[17] = m_y*R3_32*R1_33;
    expr[18] = t3_x*R1_31;  
    expr[19] = t3_y*R1_31;
    expr[20] = t3_z*R1_31;
    expr[21] = t3_x*R1_32;  
    expr[22] = t3_y*R1_32;
    expr[23] = t3_z*R1_32;
    expr[24] = t3_x*R1_33;  
    expr[25] = t3_y*R1_33;
    expr[26] = t3_z*R1_33;

    unsigned int numPoints = static_cast<unsigned int>(this->data->size());
    for(unsigned int i =0; i<numPoints; i++) { 
      lsqrRecipes::Frame &T2 = (*(this->data))[i]->T2;
      T2.getRotationMatrix(R2);
      T2.getTranslation(t2);
      u = (*(this->data))[i]->q[0];
      v = (*(this->data))[i]->q[1];
          //fx has the correct size (see vnl_least_squares_function documentation)
      fx[i] = u*R2[0][0]*expr[0] +
              u*R2[0][1]*expr[1] +
              u*R2[0][2]*expr[2] +              
              u*R2[1][0]*expr[3] +
              u*R2[1][1]*expr[4] +
              u*R2[1][2]*expr[5] +
              u*R2[2][0]*expr[6] +
              u*R2[2][1]*expr[7] +
              u*R2[2][2]*expr[8] +
              v*R2[0][0]*expr[9] +
              v*R2[0][1]*expr[10] +
              v*R2[0][2]*expr[11] +              
              v*R2[1][0]*expr[12] +
              v*R2[1][1]*expr[13] +
              v*R2[1][2]*expr[14] +
              v*R2[2][0]*expr[15] +
              v*R2[2][1]*expr[16] +
              v*R2[2][2]*expr[17] +
              R2[0][0]*expr[18] +
              R2[0][1]*expr[19] +
              R2[0][2]*expr[20] +              
              R2[1][0]*expr[21] +
              R2[1][1]*expr[22] +
              R2[1][2]*expr[23] +
              R2[2][0]*expr[24] +
              R2[2][1]*expr[25] +
              R2[2][2]*expr[26] +
              t2[0]*R1_31 +
              t2[1]*R1_32 +
              t2[2]*R1_33 +
              t1_z;
    }
}

/*****************************************************************************/
void PlanePhantomUSCalibrationParametersEstimator::SumSquaresCalibrationPointsDistanceFunction::gradf(vnl_vector<double> const& x, vnl_matrix<double>& jacobian)
{
}

} //namespace lsqrRecipes
