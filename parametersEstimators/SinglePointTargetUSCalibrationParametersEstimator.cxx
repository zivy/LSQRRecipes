#include "SinglePointTargetUSCalibrationParametersEstimator.h"
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_cross.h>
#include <vnl/algo/vnl_levenberg_marquardt.h>
#include "Epsilon.h"

namespace lsqrRecipes {

SingleUnknownPointTargetUSCalibrationParametersEstimator::SingleUnknownPointTargetUSCalibrationParametersEstimator(double delta, LeastSquaresType lsType) : 
  ParametersEstimator<DataType,double>(4) 
{
  this->deltaSquared = delta*delta;
  this->lsType = lsType;
}
/*****************************************************************************/
void SingleUnknownPointTargetUSCalibrationParametersEstimator::estimate(std::vector<DataType *> &data, 
																	                                      std::vector<double> &parameters)
{
  parameters.clear();
  if(data.size() != this->minForEstimate)
    return;
        //forward the work to the analytic least squares
  analyticLeastSquaresEstimate(data,parameters);
}
/*****************************************************************************/
void SingleUnknownPointTargetUSCalibrationParametersEstimator::estimate(std::vector<DataType> &data, 
																	                                      std::vector<double> &parameters)
{
	std::vector<DataType *> usedData;
	size_t dataSize = data.size();
	for(size_t i=0; i<dataSize; i++)
		usedData.push_back(&(data[i]));
	estimate(usedData,parameters);
}
/*****************************************************************************/
void SingleUnknownPointTargetUSCalibrationParametersEstimator::leastSquaresEstimate(std::vector<DataType *> &data, 
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
void SingleUnknownPointTargetUSCalibrationParametersEstimator::leastSquaresEstimate(std::vector<DataType> &data, 
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
bool SingleUnknownPointTargetUSCalibrationParametersEstimator::agree(std::vector<double> &parameters, 
                                                                     DataType &data)
{
  vnl_matrix<double> T2(4,4), T3(4,4);
  vnl_vector<double> q(4), qInT(4);
  double R2[3][3], t2[3];
  double errX, errY, errZ;

  data.T2.getRotationMatrix(R2);
  data.T2.getTranslation(t2);

  T2(0,0) = R2[0][0];  T2(0,1) = R2[0][1];  T2(0,2) = R2[0][2];  T2(0,3) = t2[0];
  T2(1,0) = R2[1][0];  T2(1,1) = R2[1][1];  T2(1,2) = R2[1][2];  T2(1,3) = t2[1];
  T2(2,0) = R2[2][0];  T2(2,1) = R2[2][1];  T2(2,2) = R2[2][2];  T2(2,3) = t2[2];
  T2(3,0) = 0.0;       T2(3,1) = 0.0;       T2(3,2) = 0.0;       T2(3,3) = 1;

  T3(0,0) = parameters[11];  T3(0,1) = parameters[14];  T3(0,2) = parameters[17];  T3(0,3) = parameters[3];
  T3(1,0) = parameters[12];  T3(1,1) = parameters[15];  T3(1,2) = parameters[18];  T3(1,3) = parameters[4];
  T3(2,0) = parameters[13];  T3(2,1) = parameters[16];  T3(2,2) = parameters[19];  T3(2,3) = parameters[5];
  T3(3,0) = 0.0;             T3(3,1) = 0.0;             T3(3,2) = 0.0;             T3(3,3) = 1;

  q(0) = data.q[0];
  q(1) = data.q[1];
  q(2) = 0.0;
  q(3) = 1.0;

  qInT = T2*T3*q;  

  errX = qInT(0) - parameters[0];
  errY = qInT(1) - parameters[1];
  errZ = qInT(2) - parameters[2];

  return(errX*errX+errY*errY+errZ*errZ<this->deltaSquared);
}
/*****************************************************************************/
/**
 * Solve the equation system Ax=b, where A is a 3nX12 matrix (n>=4), x is a
 * 12x1 vector and b is a 3nx1 vector, with each data element providing the
 * following three equations:
 * 
 *                                   [m_x*R3(:,1)]         
 * [u_i*R2_i  v_i*R2_i  R2_i  -I]    [m_y*R3(:,2)]    =  [-t2]
 *                                   [    t3     ]
 *                                   [    t1     ]
 *
 */
void SingleUnknownPointTargetUSCalibrationParametersEstimator::analyticLeastSquaresEstimate(std::vector< DataType *> &data, 
                                                                                             std::vector<double> &parameters)
{
	parameters.clear();
           //not enough data elements for estimate
	if(data.size()<this->minForEstimate) 
		return;

  unsigned int index, i, numPoints;

  numPoints = static_cast<unsigned int>(data.size());

  vnl_matrix<double> A(3*numPoints,12);
	vnl_vector<double> x(12);
	vnl_vector<double> b(3*numPoints); 
  double R2[3][3], t2x, t2y, t2z, ui, vi;

  for(i=0; i<numPoints; i++) {
           //transformation from US reference frame to tracker
    lsqrRecipes::Frame &frm = data[i]->T2;               
    frm.getRotationMatrix(R2);
    frm.getTranslation(t2x,t2y,t2z);
           //2D pixel coordinates
    ui = data[i]->q[0];
    vi = data[i]->q[1];
    index = 3*i; 
            //first row for current data element
    A(index,0) = R2[0][0]*ui; 
    A(index,1) = R2[0][1]*ui;
    A(index,2) = R2[0][2]*ui;
    A(index,3) = R2[0][0]*vi; 
    A(index,4) = R2[0][1]*vi;
    A(index,5) = R2[0][2]*vi;
    A(index,6) = R2[0][0]; 
    A(index,7) = R2[0][1];
    A(index,8) = R2[0][2];
    A(index,9) = -1.0; 
    A(index,10) = 0.0;
    A(index,11) = 0.0;
    b(index) = -t2x;
    index++;
            //second row for current data element
    A(index,0) = R2[1][0]*ui; 
    A(index,1) = R2[1][1]*ui;
    A(index,2) = R2[1][2]*ui;
    A(index,3) = R2[1][0]*vi; 
    A(index,4) = R2[1][1]*vi;
    A(index,5) = R2[1][2]*vi;
    A(index,6) = R2[1][0]; 
    A(index,7) = R2[1][1];
    A(index,8) = R2[1][2];
    A(index,9) = 0.0; 
    A(index,10) = -1.0;
    A(index,11) = 0.0;
    b(index) = -t2y;
    index++;
            //third row for current data element
    A(index,0) = R2[2][0]*ui; 
    A(index,1) = R2[2][1]*ui;
    A(index,2) = R2[2][2]*ui;
    A(index,3) = R2[2][0]*vi; 
    A(index,4) = R2[2][1]*vi;
    A(index,5) = R2[2][2]*vi;
    A(index,6) = R2[2][0]; 
    A(index,7) = R2[2][1];
    A(index,8) = R2[2][2];
    A(index,9) = 0.0; 
    A(index,10) = 0.0;
    A(index,11) = -1.0;
    b(index) = -t2z;
  }

  vnl_matrix_inverse<double> Ainv(A);
               //explicitly zero out small singular values
               //this is ugly as it exposes that the inverse is computed via SVD
               //the magic number I use as threshold is FLT_EPSILON (see float.h)
  double onlyTranslationEPS = 1.192092896e-07;
  Ainv.zero_out_absolute(onlyTranslationEPS);

	if(Ainv.rank()<12) //points do not yield a solution
	  return;
	x = Ainv * b;
  
        //get the scale factors and rotation angles
  double omega_z, omega_y, omega_x, m_x, m_y;
  vnl_vector<double> r1(3), r2(3), r3(3);
  vnl_matrix<double> R3(3,3);
       //get the x scale factor and the first column of the rotation matrix
  r1(0) = x(0);    
  r1(1) = x(1);    
  r1(2) = x(2);
  m_x = r1.two_norm();
  r1.normalize();
       //get the y scale factor and the second column of the rotation matrix
  r2(0) = x(3);    
  r2(1) = x(4);    
  r2(2) = x(5);
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
  double smallAngle = 0.008726535498373935;//0.5 degrees
  double halfPI = 1.5707963267948966192313216916398;

           //two options for omega_y depending on choice of +-1 for sqrt
           //we arbitrarily choose +1
  omega_y = atan2(-R3(2,0), sqrt(R3(0,0)*R3(0,0) + R3(1,0)*R3(1,0)));
             //omega_z and omega_x depend on the omega_y value
	if(fabs(omega_y - halfPI) > smallAngle && fabs(omega_y + halfPI) > smallAngle) {
    double cy = cos(omega_y);
    omega_z = atan2(R3(1,0)/cy, R3(0,0)/cy);
    omega_x = atan2(R3(2,1)/cy, R3(2,2)/cy);
  }     //gimbal lock, omega_y is approximatly plus/minus half PI, arbitrarily 
        //set omega_z to 0 and compute omega_x
  else {          
    omega_z = 0;
    omega_x = atan2(R3(0,1), R3(1,1));
  }

       //set the parameters    
  parameters.push_back(x(9));  //t1_x 
  parameters.push_back(x(10)); //t1_y
  parameters.push_back(x(11)); //t1_z
  parameters.push_back(x(6));  //t3_x 
  parameters.push_back(x(7)); //t3_y
  parameters.push_back(x(8)); //t3_z
  parameters.push_back(omega_z);
  parameters.push_back(omega_y);
  parameters.push_back(omega_x);
  parameters.push_back(m_x);
  parameters.push_back(m_y);
  parameters.push_back(m_x*R3(0,0));
  parameters.push_back(m_x*R3(1,0));
  parameters.push_back(m_x*R3(2,0));
  parameters.push_back(m_y*R3(0,1));
  parameters.push_back(m_y*R3(1,1));
  parameters.push_back(m_y*R3(2,1));
  parameters.push_back(R3(0,2));
  parameters.push_back(R3(1,2));
  parameters.push_back(R3(2,2));
}
/*****************************************************************************/
void SingleUnknownPointTargetUSCalibrationParametersEstimator::iterativeLeastSquaresEstimate(std::vector< DataType *> &data, 
																	                                                           std::vector<double> &initialParameters, 
                                                                                             std::vector<double> &finalParameters)
{
  unsigned int i, numParameters = 11;
	vnl_vector<double> parameters(numParameters);
	
	for(i=0; i<numParameters; i++)
	  parameters[i] = initialParameters[i];

	SingleUnknownPointTargetUSCalibrationParametersEstimator::SumSquaresCalibrationPointsDistanceFunction 
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
    double cz,sz,cy,sy,cx,sx;
    cz = cos(finalParameters[6]);
    sz = sin(finalParameters[6]);
    cy = cos(finalParameters[7]);
    sy = sin(finalParameters[7]);
    cx = cos(finalParameters[8]);
    sx = sin(finalParameters[8]);
         //m_x*R3_11
    finalParameters.push_back(finalParameters[9]*cz*cy);
         //m_x*R3_21
    finalParameters.push_back(finalParameters[9]*sz*cy);
         //m_x*R3_31
    finalParameters.push_back(-finalParameters[9]*sy);
         //m_y*R3_12
    finalParameters.push_back(finalParameters[10]*(cz*sy*sx-sz*cx));
         //m_y*R3_22
    finalParameters.push_back(finalParameters[10]*(sz*sy*sx+cz*cx));
         //m_y*R3_32
    finalParameters.push_back(finalParameters[10]*cy*sx);
         //R3_13
    finalParameters.push_back(cz*sy*cx+sz*sx);
         //R3_23
    finalParameters.push_back(sz*sy*cx-cz*sx);
         //R3_33
    finalParameters.push_back(cy*cx);
  }
}
/*****************************************************************************/
void SingleUnknownPointTargetUSCalibrationParametersEstimator::getDistanceStatistics(const std::vector<double> &parameters, 
                                                                                     const std::vector<DataType> &data, 
                                                                                     std::vector<double> &distances,
                                                                                     double &min, double &max, double &mean)
{
  vnl_matrix<double> T2(4,4), T3(4,4);
  vnl_vector<double> q(4), qInT(4);
  double R2[3][3], t2[3];
  double errX, errY, errZ, dist;
  size_t i, n;

            //number of parameters doesn't match expected parameter vector
            //[t_1x, t_1y, t_1z, t_3x, t_3y, t_3z, omega_z, omega_y, omega_x, m_x, m_y, m_x*R3(:,1), m_y*R3(:,2), R3(:,3)]
  if(parameters.size()<20)
    throw std::exception();

  T3(0,0) = parameters[11];  T3(0,1) = parameters[14];  T3(0,2) = parameters[17];  T3(0,3) = parameters[3];
  T3(1,0) = parameters[12];  T3(1,1) = parameters[15];  T3(1,2) = parameters[18];  T3(1,3) = parameters[4];
  T3(2,0) = parameters[13];  T3(2,1) = parameters[16];  T3(2,2) = parameters[19];  T3(2,3) = parameters[5];
  T3(3,0) = 0.0;             T3(3,1) = 0.0;             T3(3,2) = 0.0;             T3(3,3) = 1;


          //use first point for initialization
  data[0].T2.getRotationMatrix(R2);
  data[0].T2.getTranslation(t2);
  T2(0,0) = R2[0][0];  T2(0,1) = R2[0][1];  T2(0,2) = R2[0][2];  T2(0,3) = t2[0];
  T2(1,0) = R2[1][0];  T2(1,1) = R2[1][1];  T2(1,2) = R2[1][2];  T2(1,3) = t2[1];
  T2(2,0) = R2[2][0];  T2(2,1) = R2[2][1];  T2(2,2) = R2[2][2];  T2(2,3) = t2[2];
  T2(3,0) = 0.0;       T2(3,1) = 0.0;       T2(3,2) = 0.0;       T2(3,3) = 1;

  q(0) = data[0].q[0];
  q(1) = data[0].q[1];
  q(2) = 0.0;
  q(3) = 1.0;

  qInT = T2*T3*q;  
  errX = qInT(0) - parameters[0];
  errY = qInT(1) - parameters[1];
  errZ = qInT(2) - parameters[2];
  distances.push_back(sqrt(errX*errX+errY*errY+errZ*errZ));
  min = max = mean = distances[0];
             
  n = data.size();
            //go over the rest of the points
  for(i=1; i<n; i++) {
    data[i].T2.getRotationMatrix(R2);
    data[i].T2.getTranslation(t2);
    T2(0,0) = R2[0][0];  T2(0,1) = R2[0][1];  T2(0,2) = R2[0][2];  T2(0,3) = t2[0];
    T2(1,0) = R2[1][0];  T2(1,1) = R2[1][1];  T2(1,2) = R2[1][2];  T2(1,3) = t2[1];
    T2(2,0) = R2[2][0];  T2(2,1) = R2[2][1];  T2(2,2) = R2[2][2];  T2(2,3) = t2[2];
    T2(3,0) = 0.0;       T2(3,1) = 0.0;       T2(3,2) = 0.0;       T2(3,3) = 1;

    q(0) = data[i].q[0];
    q(1) = data[i].q[1];
    q(2) = 0.0;
    q(3) = 1.0;

    qInT = T2*T3*q;  
    errX = qInT(0) - parameters[0];
    errY = qInT(1) - parameters[1];
    errZ = qInT(2) - parameters[2];
    dist = sqrt(errX*errX+errY*errY+errZ*errZ);
    distances.push_back(dist);
    mean+=dist;
    if(dist>max)
      max = dist;
    else if(dist<min)
      min = dist;
  }
  mean/= n;
}
/*****************************************************************************/
SingleUnknownPointTargetUSCalibrationParametersEstimator::SumSquaresCalibrationPointsDistanceFunction::SumSquaresCalibrationPointsDistanceFunction(std::vector<DataType *> *data) :
vnl_least_squares_function(11, static_cast<unsigned int>(data->size()), vnl_least_squares_function::use_gradient)	                                                  
{
	this->data = data;
}

/*****************************************************************************/
void SingleUnknownPointTargetUSCalibrationParametersEstimator::SumSquaresCalibrationPointsDistanceFunction::setData(std::vector<DataType *> *data)
{
	this->data = data;																											
}
/*****************************************************************************/
void SingleUnknownPointTargetUSCalibrationParametersEstimator::SumSquaresCalibrationPointsDistanceFunction::f(vnl_vector<double> const &x, vnl_vector<double> &fx)
{
  double t_1x, t_1y, t_1z, t_3x, t_3y, t_3z, sz, cz, sy, cy, sx, cx, m_x, m_y;
  double R3_11, R3_21, R3_31, R3_12, R3_22, R3_32;

  t_1x = x(0);
  t_1y = x(1);
  t_1z = x(2);
  t_3x = x(3);
  t_3y = x(4); 
  t_3z = x(5);
  sz = sin(x(6));
  cz = cos(x(6));
  sy = sin(x(7));
  cy = cos(x(7));
  sx = sin(x(8));
  cx = cos(x(8));
  m_x = x(9);
  m_y = x(10);

  R3_11 = cz*cy;
  R3_21 = sz*cy;
  R3_31 = -sy;
  R3_12 = cz*sy*sx - sz*cx;
  R3_22 = sz*sy*sx+cz*cx;
  R3_32 = cy*sx;
 
  unsigned int numPoints = static_cast<unsigned int>(this->data->size());
  double R2[3][3], t2[3];
  double A_11, A_12, A_13, A_14, A_15, A_16, A_17, A_18, A_19;
  double A_21, A_22, A_23, A_24, A_25, A_26, A_27, A_28, A_29;
  double A_31, A_32, A_33, A_34, A_35, A_36, A_37, A_38, A_39;
  double b_1, b_2, b_3, expr1, expr2, expr3;

  for(unsigned int i =0; i<numPoints; i++) { 
    lsqrRecipes::Frame &T2 = (*(this->data))[i]->T2;
    T2.getRotationMatrix(R2);
    T2.getTranslation(t2);
    lsqrRecipes::Point2D &q = (*(this->data))[i]->q;
             //u_i*R_2i
    A_11 = q[0]*R2[0][0];   A_12 = q[0]*R2[0][1];  A_13 = q[0]*R2[0][2];
    A_21 = q[0]*R2[1][0];   A_22 = q[0]*R2[1][1];  A_23 = q[0]*R2[1][2];
    A_31 = q[0]*R2[2][0];   A_32 = q[0]*R2[2][1];  A_33 = q[0]*R2[2][2];
             //v_i*R_2i
    A_14 = q[1]*R2[0][0];   A_15 = q[1]*R2[0][1];  A_16 = q[1]*R2[0][2];
    A_24 = q[1]*R2[1][0];   A_25 = q[1]*R2[1][1];  A_26 = q[1]*R2[1][2];
    A_34 = q[1]*R2[2][0];   A_35 = q[1]*R2[2][1];  A_36 = q[1]*R2[2][2];
               //R_2i
    A_17 = R2[0][0];   A_18 = R2[0][1];  A_19 = R2[0][2];
    A_27 = R2[1][0];   A_28 = R2[1][1];  A_29 = R2[1][2];
    A_37 = R2[2][0];   A_38 = R2[2][1];  A_39 = R2[2][2];

    b_1 = -t2[0];
    b_2 = -t2[1];
    b_3 = -t2[2];

    expr1 = A_11*m_x*R3_11 +
            A_12*m_x*R3_21 +
            A_13*m_x*R3_31 + 
            A_14*m_y*R3_12 + 
            A_15*m_y*R3_22 + 
            A_16*m_y*R3_32 + 
            A_17*t_3x +
            A_18*t_3y +
            A_19*t_3z -
            t_1x - 
            b_1;
        
    expr2 = A_21*m_x*R3_11 + 
            A_22*m_x*R3_21 + 
            A_23*m_x*R3_31 + 
            A_24*m_y*R3_12 +
            A_25*m_y*R3_22 +
            A_26*m_y*R3_32 +
            A_27*t_3x +
            A_28*t_3y +
            A_29*t_3z -
            t_1y -
            b_2;

    expr3 = A_31*m_x*R3_11 +
            A_32*m_x*R3_21 +
            A_33*m_x*R3_31 +
            A_34*m_y*R3_12 +
            A_35*m_y*R3_22 +
            A_36*m_y*R3_32 +
            A_37*t_3x +
            A_38*t_3y +
            A_39*t_3z -
            t_1z -
            b_3;
          //fx has the correct size (see vnl_least_squares_function documentation)
    fx[i] = sqrt(expr1*expr1+expr2*expr2+expr3*expr3);
  }
}

/*****************************************************************************/
void SingleUnknownPointTargetUSCalibrationParametersEstimator::SumSquaresCalibrationPointsDistanceFunction::gradf(vnl_vector<double> const& x, vnl_matrix<double>& jacobian)
{
  double t_1x, t_1y, t_1z, t_3x, t_3y, t_3z, sz, cz, sy, cy, sx, cx, m_x, m_y;
  double R3_11, R3_21, R3_31, R3_12, R3_22, R3_32;

  t_1x = x(0);
  t_1y = x(1);
  t_1z = x(2);
  t_3x = x(3);
  t_3y = x(4); 
  t_3z = x(5);
  sz = sin(x(6));
  cz = cos(x(6));
  sy = sin(x(7));
  cy = cos(x(7));
  sx = sin(x(8));
  cx = cos(x(8));
  m_x = x(9);
  m_y = x(10);

  R3_11 = cz*cy;
  R3_21 = sz*cy;
  R3_31 = -sy;
  R3_12 = cz*sy*sx - sz*cx;
  R3_22 = sz*sy*sx+cz*cx;
  R3_32 = cy*sx;

  unsigned int numPoints = static_cast<unsigned int>(this->data->size());
  double R2[3][3], t2[3];
  double A_11, A_12, A_13, A_14, A_15, A_16, A_17, A_18, A_19;
  double A_21, A_22, A_23, A_24, A_25, A_26, A_27, A_28, A_29;
  double A_31, A_32, A_33, A_34, A_35, A_36, A_37, A_38, A_39;
  double b_1, b_2, b_3, expr1, expr2, expr3, delta_i;
  double v1, v2, v3, v4, v5, v6;

  for(unsigned int i =0; i<numPoints; i++) { 
    lsqrRecipes::Frame &T2 = (*(this->data))[i]->T2;
    T2.getRotationMatrix(R2);
    T2.getTranslation(t2);
    lsqrRecipes::Point2D &q = (*(this->data))[i]->q;
             //u_i*R_2i
    A_11 = q[0]*R2[0][0];   A_12 = q[0]*R2[0][1];  A_13 = q[0]*R2[0][2];
    A_21 = q[0]*R2[1][0];   A_22 = q[0]*R2[1][1];  A_23 = q[0]*R2[1][2];
    A_31 = q[0]*R2[2][0];   A_32 = q[0]*R2[2][1];  A_33 = q[0]*R2[2][2];
             //v_i*R_2i
    A_14 = q[1]*R2[0][0];   A_15 = q[1]*R2[0][1];  A_16 = q[1]*R2[0][2];
    A_24 = q[1]*R2[1][0];   A_25 = q[1]*R2[1][1];  A_26 = q[1]*R2[1][2];
    A_34 = q[1]*R2[2][0];   A_35 = q[1]*R2[2][1];  A_36 = q[1]*R2[2][2];
               //R_2i
    A_17 = R2[0][0];   A_18 = R2[0][1];  A_19 = R2[0][2];
    A_27 = R2[1][0];   A_28 = R2[1][1];  A_29 = R2[1][2];
    A_37 = R2[2][0];   A_38 = R2[2][1];  A_39 = R2[2][2];

    b_1 = -t2[0];
    b_2 = -t2[1];
    b_3 = -t2[2];

    expr1 = A_11*m_x*R3_11 +
            A_12*m_x*R3_21 +
            A_13*m_x*R3_31 + 
            A_14*m_y*R3_12 + 
            A_15*m_y*R3_22 + 
            A_16*m_y*R3_32 + 
            A_17*t_3x +
            A_18*t_3y +
            A_19*t_3z -
            t_1x - 
            b_1;
        
    expr2 = A_21*m_x*R3_11 + 
            A_22*m_x*R3_21 + 
            A_23*m_x*R3_31 + 
            A_24*m_y*R3_12 +
            A_25*m_y*R3_22 +
            A_26*m_y*R3_32 +
            A_27*t_3x +
            A_28*t_3y +
            A_29*t_3z -
            t_1y -
            b_2;

    expr3 = A_31*m_x*R3_11 +
            A_32*m_x*R3_21 +
            A_33*m_x*R3_31 +
            A_34*m_y*R3_12 +
            A_35*m_y*R3_22 +
            A_36*m_y*R3_32 +
            A_37*t_3x +
            A_38*t_3y +
            A_39*t_3z -
            t_1z -
            b_3;

    delta_i = sqrt(expr1*expr1+expr2*expr2+expr3*expr3);     
   
             //partial derivative t_1x
    jacobian(i,0) = -expr1/delta_i;
             //partial derivative t_1y   
    jacobian(i,1) = -expr2/delta_i;
             //partial derivative t_1z   
    jacobian(i,2) = -expr3/delta_i;
             //partial derivative t_3x   
    jacobian(i,3) = (A_17*expr1+A_27*expr2+A_37*expr3)/delta_i;
                //partial derivative t_3y   
    jacobian(i,4) = (A_18*expr1+A_28*expr2+A_38*expr3)/delta_i;
                //partial derivative t_3z   
    jacobian(i,5) = (A_19*expr1+A_29*expr2+A_39*expr3)/delta_i;   
         //partial derivative omega_z
    v1 = -m_x*sz*cy;
    v2 = m_x*cz*cy;
    v3 = -m_y*(sz*sy*sx+cz*cx);
    v4 = m_y*(cz*sy*sx-sz*cx);   
    jacobian(i,6) = ((A_11*v1+A_12*v2+A_14*v3+A_15*v4)*expr1 +
                     (A_21*v1+A_22*v2+A_24*v3+A_25*v4)*expr2 +
                     (A_31*v1+A_32*v2+A_34*v3+A_35*v4)*expr3)/
                     delta_i;
          //partial derivative omega_y
    v1 = -m_x*sy*cz;
    v2 = -m_x*sy*sz;
    v3 = -m_x*cy;
    v4 = m_y*sx*cy*cz;
    v5 = m_y*sx*cy*sz;
    v6 = -m_y*sx*sy;
    jacobian(i,7) = ((A_11*v1+A_12*v2+A_13*v3+A_14*v4+A_15*v5+A_16*v6)*expr1 +
                     (A_21*v1+A_22*v2+A_23*v3+A_24*v4+A_25*v5+A_26*v6)*expr2 +
                     (A_31*v1+A_32*v2+A_33*v3+A_34*v4+A_35*v5+A_36*v6)*expr3)/
                     delta_i;                  
          //partial derivative omega_x
    v1 = m_y*(cz*sy*cx+sz*sx);
    v2 = m_y*(sz*sy*cx-cz*sx);
    v3 = m_y*cy*cx;
    jacobian(i,8) = ((A_14*v1+A_15*v2+A_16*v3)*expr1 +
                     (A_24*v1+A_25*v2+A_26*v3)*expr2 +
                     (A_34*v1+A_35*v2+A_36*v3)*expr3)/
                     delta_i;                  
             //partial derivative m_x
    jacobian(i,9) = ((A_11*R3_11+A_12*R3_21+A_13*R3_31)*expr1 +
                      (A_21*R3_11+A_22*R3_21+A_23*R3_31)*expr2 +
                      (A_31*R3_11+A_32*R3_21+A_33*R3_31)*expr3)/
                      delta_i;
             //partial derivative m_y                   
    jacobian(i,10) = ((A_14*R3_12+A_15*R3_22+A_16*R3_32)*expr1 +
                      (A_24*R3_12+A_25*R3_22+A_26*R3_32)*expr2 +
                      (A_34*R3_12+A_35*R3_22+A_36*R3_32)*expr3)/
                      delta_i;                                      
  }
}


/*****************************************************************************/
/*****************************************************************************/

CalibratedPointerTargetUSCalibrationParametersEstimator::CalibratedPointerTargetUSCalibrationParametersEstimator(double delta, LeastSquaresType lsType) : 
  ParametersEstimator<DataType,double>(3) 
{
  this->deltaSquared = delta*delta;
  this->lsType = lsType;
}
/*****************************************************************************/
void CalibratedPointerTargetUSCalibrationParametersEstimator::estimate(std::vector<DataType *> &data, 
																	                                      std::vector<double> &parameters)
{
  parameters.clear();
  if(data.size() != this->minForEstimate)
    return;
        //forward the work to the analytic least squares
  analyticLeastSquaresEstimate(data,parameters);
}
/*****************************************************************************/
void CalibratedPointerTargetUSCalibrationParametersEstimator::estimate(std::vector<DataType> &data, 
																	                                      std::vector<double> &parameters)
{
	std::vector<DataType *> usedData;
	size_t dataSize = data.size();
	for(size_t i=0; i<dataSize; i++)
		usedData.push_back(&(data[i]));
	estimate(usedData,parameters);
}
/*****************************************************************************/
void CalibratedPointerTargetUSCalibrationParametersEstimator::leastSquaresEstimate(std::vector<DataType *> &data, 
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
void CalibratedPointerTargetUSCalibrationParametersEstimator::leastSquaresEstimate(std::vector<DataType> &data, 
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
bool CalibratedPointerTargetUSCalibrationParametersEstimator::agree(std::vector<double> &parameters, 
                                                                     DataType &data)
{
  vnl_matrix<double> T2(4,4), T3(4,4);
  vnl_vector<double> q(4), qInT(4);
  double R2[3][3], t2[3];
  double errX, errY, errZ;

  data.T2.getRotationMatrix(R2);
  data.T2.getTranslation(t2);

  T2(0,0) = R2[0][0];  T2(0,1) = R2[0][1];  T2(0,2) = R2[0][2];  T2(0,3) = t2[0];
  T2(1,0) = R2[1][0];  T2(1,1) = R2[1][1];  T2(1,2) = R2[1][2];  T2(1,3) = t2[1];
  T2(2,0) = R2[2][0];  T2(2,1) = R2[2][1];  T2(2,2) = R2[2][2];  T2(2,3) = t2[2];
  T2(3,0) = 0.0;       T2(3,1) = 0.0;       T2(3,2) = 0.0;       T2(3,3) = 1;

  T3(0,0) = parameters[8];  T3(0,1) = parameters[11];  T3(0,2) = parameters[14];  T3(0,3) = parameters[0];
  T3(1,0) = parameters[9];  T3(1,1) = parameters[12];  T3(1,2) = parameters[15];  T3(1,3) = parameters[1];
  T3(2,0) = parameters[10]; T3(2,1) = parameters[13];  T3(2,2) = parameters[16];  T3(2,3) = parameters[2];
  T3(3,0) = 0.0;            T3(3,1) = 0.0;             T3(3,2) = 0.0;             T3(3,3) = 1;

  q(0) = data.q[0];
  q(1) = data.q[1];
  q(2) = 0.0;
  q(3) = 1.0;

  qInT = T2*T3*q;  

  errX = qInT(0) - data.p[0];
  errY = qInT(1) - data.p[1];
  errZ = qInT(2) - data.p[2];

  return(errX*errX+errY*errY+errZ*errZ<this->deltaSquared);
}
/*****************************************************************************/
/**
 * Solve the equation system Ax=b, where A is a 3nX9 matrix (n>=3), x is a
 * 9x1 vector and b is a 3nx1 vector, with each data element providing the
 * following three equations:
 * 
 *                               [m_x*R3(:,1)]         
 * [u_i*R2_i  v_i*R2_i  R2_i]    [m_y*R3(:,2)]    =  [p_i-t2_i]
 *                               [    t3     ]
 *
 */
void CalibratedPointerTargetUSCalibrationParametersEstimator::analyticLeastSquaresEstimate(std::vector< DataType *> &data, 
                                                                                           std::vector<double> &parameters)
{
	parameters.clear();
           //not enough data elements for estimate
	if(data.size()<this->minForEstimate) 
		return;

  unsigned int index, i, numPoints;

  numPoints = static_cast<unsigned int>(data.size());
  vnl_matrix<double> A(3*numPoints,9);
	vnl_vector<double> x(9);
	vnl_vector<double> b(3*numPoints); 
  double R2[3][3], t2x, t2y, t2z, ui, vi, pi_x, pi_y, pi_z;
  
  for(i=0; i<numPoints; i++) {
           //transformation from US reference frame to tracker
    lsqrRecipes::Frame &frm = data[i]->T2;               
    frm.getRotationMatrix(R2);
    frm.getTranslation(t2x,t2y,t2z);
           //2D pixel coordinates
    ui = data[i]->q[0];
    vi = data[i]->q[1];
          //3D pointer tip coordinates
    pi_x = data[i]->p[0];
    pi_y = data[i]->p[1];
    pi_z = data[i]->p[2];
    index = 3*i; 
            //first row for current data element
    A(index,0) = R2[0][0]*ui; 
    A(index,1) = R2[0][1]*ui;
    A(index,2) = R2[0][2]*ui;
    A(index,3) = R2[0][0]*vi; 
    A(index,4) = R2[0][1]*vi;
    A(index,5) = R2[0][2]*vi;
    A(index,6) = R2[0][0]; 
    A(index,7) = R2[0][1];
    A(index,8) = R2[0][2];
    b(index) = pi_x-t2x;
    index++;
            //second row for current data element
    A(index,0) = R2[1][0]*ui; 
    A(index,1) = R2[1][1]*ui;
    A(index,2) = R2[1][2]*ui;
    A(index,3) = R2[1][0]*vi; 
    A(index,4) = R2[1][1]*vi;
    A(index,5) = R2[1][2]*vi;
    A(index,6) = R2[1][0]; 
    A(index,7) = R2[1][1];
    A(index,8) = R2[1][2];
    b(index) = pi_y-t2y;
    index++;
            //third row for current data element
    A(index,0) = R2[2][0]*ui; 
    A(index,1) = R2[2][1]*ui;
    A(index,2) = R2[2][2]*ui;
    A(index,3) = R2[2][0]*vi; 
    A(index,4) = R2[2][1]*vi;
    A(index,5) = R2[2][2]*vi;
    A(index,6) = R2[2][0]; 
    A(index,7) = R2[2][1];
    A(index,8) = R2[2][2];
    b(index) = pi_z-t2z;
  }

  vnl_matrix_inverse<double> Ainv(A);
               //explicitly zero out small singular values
               //this is ugly as it exposes that the inverse is computed via SVD
               //the magic number I use as threshold is FLT_EPSILON (see float.h)
  double singularValueEPS = 1.192092896e-07;
  Ainv.zero_out_absolute(singularValueEPS);

	if(Ainv.rank()<9) //points do not yield a solution
	  return;
	x = Ainv * b;
  
        //get the scale factors and rotation angles
  double omega_z, omega_y, omega_x, m_x, m_y;
  vnl_vector<double> r1(3), r2(3), r3(3);
  vnl_matrix<double> R3(3,3);
       //get the x scale factor and the first column of the rotation matrix
  r1(0) = x(0);    
  r1(1) = x(1);    
  r1(2) = x(2);
  m_x = r1.two_norm();
  r1.normalize();
       //get the y scale factor and the second column of the rotation matrix
  r2(0) = x(3);    
  r2(1) = x(4);    
  r2(2) = x(5);
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
  double smallAngle = 0.008726535498373935;//0.5 degrees
  double halfPI = 1.5707963267948966192313216916398;

           //two options for omega_y depending on choice of +-1 for sqrt
           //we arbitrarily choose +1
  omega_y = atan2(-R3(2,0), sqrt(R3(0,0)*R3(0,0) + R3(1,0)*R3(1,0)));
             //omega_z and omega_x depend on the omega_y value
	if(fabs(omega_y - halfPI) > smallAngle && fabs(omega_y + halfPI) > smallAngle) {
    double cy = cos(omega_y);
    omega_z = atan2(R3(1,0)/cy, R3(0,0)/cy);
    omega_x = atan2(R3(2,1)/cy, R3(2,2)/cy);
  }     //gimbal lock, omega_y is approximatly plus/minus half PI, arbitrarily 
        //set omega_z to 0 and compute omega_x
  else {          
    omega_z = 0;
    omega_x = atan2(R3(0,1), R3(1,1));
  }

       //set the parameters    
  parameters.push_back(x(6));  //t3_x 
  parameters.push_back(x(7)); //t3_y
  parameters.push_back(x(8)); //t3_z
  parameters.push_back(omega_z);
  parameters.push_back(omega_y);
  parameters.push_back(omega_x);
  parameters.push_back(m_x);
  parameters.push_back(m_y);
  parameters.push_back(m_x*R3(0,0));
  parameters.push_back(m_x*R3(1,0));
  parameters.push_back(m_x*R3(2,0));
  parameters.push_back(m_y*R3(0,1));
  parameters.push_back(m_y*R3(1,1));
  parameters.push_back(m_y*R3(2,1));
  parameters.push_back(R3(0,2));
  parameters.push_back(R3(1,2));
  parameters.push_back(R3(2,2));
}
/*****************************************************************************/
void CalibratedPointerTargetUSCalibrationParametersEstimator::iterativeLeastSquaresEstimate(std::vector< DataType *> &data, 
												  				                                                          std::vector<double> &initialParameters, 
                                                                                            std::vector<double> &finalParameters)
{
  unsigned int i, numParameters = 8;
	vnl_vector<double> parameters(numParameters);
	
	for(i=0; i<numParameters; i++)
	  parameters[i] = initialParameters[i];

	CalibratedPointerTargetUSCalibrationParametersEstimator::SumSquaresCalibrationPointsDistanceFunction 
    optimizedFunction(&data);		
  vnl_levenberg_marquardt lmOptimization(optimizedFunction);

            //set all the optimization tolerances (see vnl_nonlinear_minimizer)
	double gradTolerance = 10e-8;
	double parametersChangeTolerance = 10e-8;
  double functionChangeTolerance = 10e-8;
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
    double cz,sz,cy,sy,cx,sx;
    cz = cos(finalParameters[3]);
    sz = sin(finalParameters[3]);
    cy = cos(finalParameters[4]);
    sy = sin(finalParameters[4]);
    cx = cos(finalParameters[5]);
    sx = sin(finalParameters[5]);
         //m_x*R3_11
    finalParameters.push_back(finalParameters[6]*cz*cy);
         //m_x*R3_21
    finalParameters.push_back(finalParameters[6]*sz*cy);
         //m_x*R3_31
    finalParameters.push_back(-finalParameters[6]*sy);
         //m_y*R3_12
    finalParameters.push_back(finalParameters[7]*(cz*sy*sx-sz*cx));
         //m_y*R3_22
    finalParameters.push_back(finalParameters[7]*(sz*sy*sx+cz*cx));
         //m_y*R3_32
    finalParameters.push_back(finalParameters[7]*cy*sx);
         //R3_13
    finalParameters.push_back(cz*sy*cx+sz*sx);
         //R3_23
    finalParameters.push_back(sz*sy*cx-cz*sx);
         //R3_33
    finalParameters.push_back(cy*cx);
  }
}
/*****************************************************************************/
void CalibratedPointerTargetUSCalibrationParametersEstimator::getDistanceStatistics(const std::vector<double> &parameters, 
                                                                                     const std::vector<DataType> &data, 
                                                                                     std::vector<double> &distances,
                                                                                     double &min, double &max, double &mean)
{
  vnl_matrix<double> T2(4,4), T3(4,4);
  vnl_vector<double> q(4), qInT(4);
  double R2[3][3], t2[3];
  double errX, errY, errZ, dist;
  size_t i, n;

            //number of parameters doesn't match expected parameter vector
            //[t_3x, t_3y, t_3z, omega_z, omega_y, omega_x, m_x, m_y, m_x*R3(:,1), m_y*R3(:,2), R3(:,3)]
  if(parameters.size()<17)
    throw std::exception();

  T3(0,0) = parameters[8];  T3(0,1) = parameters[11];  T3(0,2) = parameters[14];  T3(0,3) = parameters[0];
  T3(1,0) = parameters[9];  T3(1,1) = parameters[12];  T3(1,2) = parameters[15];  T3(1,3) = parameters[1];
  T3(2,0) = parameters[10]; T3(2,1) = parameters[13];  T3(2,2) = parameters[16];  T3(2,3) = parameters[2];
  T3(3,0) = 0.0;            T3(3,1) = 0.0;             T3(3,2) = 0.0;             T3(3,3) = 1;


          //use first point for initialization
  data[0].T2.getRotationMatrix(R2);
  data[0].T2.getTranslation(t2);
  T2(0,0) = R2[0][0];  T2(0,1) = R2[0][1];  T2(0,2) = R2[0][2];  T2(0,3) = t2[0];
  T2(1,0) = R2[1][0];  T2(1,1) = R2[1][1];  T2(1,2) = R2[1][2];  T2(1,3) = t2[1];
  T2(2,0) = R2[2][0];  T2(2,1) = R2[2][1];  T2(2,2) = R2[2][2];  T2(2,3) = t2[2];
  T2(3,0) = 0.0;       T2(3,1) = 0.0;       T2(3,2) = 0.0;       T2(3,3) = 1;

  q(0) = data[0].q[0];
  q(1) = data[0].q[1];
  q(2) = 0.0;
  q(3) = 1.0;

  qInT = T2*T3*q;  
  errX = qInT(0) - data[0].p[0];
  errY = qInT(1) - data[0].p[1];
  errZ = qInT(2) - data[0].p[2];
  distances.push_back(sqrt(errX*errX+errY*errY+errZ*errZ));
  min = max = mean = distances[0];
             
  n = data.size();
            //go over the rest of the points
  for(i=1; i<n; i++) {
    data[i].T2.getRotationMatrix(R2);
    data[i].T2.getTranslation(t2);
    T2(0,0) = R2[0][0];  T2(0,1) = R2[0][1];  T2(0,2) = R2[0][2];  T2(0,3) = t2[0];
    T2(1,0) = R2[1][0];  T2(1,1) = R2[1][1];  T2(1,2) = R2[1][2];  T2(1,3) = t2[1];
    T2(2,0) = R2[2][0];  T2(2,1) = R2[2][1];  T2(2,2) = R2[2][2];  T2(2,3) = t2[2];
    T2(3,0) = 0.0;       T2(3,1) = 0.0;       T2(3,2) = 0.0;       T2(3,3) = 1;

    q(0) = data[i].q[0];
    q(1) = data[i].q[1];
    q(2) = 0.0;
    q(3) = 1.0;

    qInT = T2*T3*q;  
    errX = qInT(0) - data[i].p[0];
    errY = qInT(1) - data[i].p[1];
    errZ = qInT(2) - data[i].p[2];
    dist = sqrt(errX*errX+errY*errY+errZ*errZ);
    distances.push_back(dist);
    mean+=dist;
    if(dist>max)
      max = dist;
    else if(dist<min)
      min = dist;
  }
  mean/= n;
}
/*****************************************************************************/
CalibratedPointerTargetUSCalibrationParametersEstimator::SumSquaresCalibrationPointsDistanceFunction::SumSquaresCalibrationPointsDistanceFunction(std::vector<DataType *> *data) :
vnl_least_squares_function(8, static_cast<unsigned int>(data->size()), vnl_least_squares_function::use_gradient)	                                                  
{
	this->data = data;
}

/*****************************************************************************/
void CalibratedPointerTargetUSCalibrationParametersEstimator::SumSquaresCalibrationPointsDistanceFunction::setData(std::vector<DataType *> *data)
{
	this->data = data;																											
}
/*****************************************************************************/
void CalibratedPointerTargetUSCalibrationParametersEstimator::SumSquaresCalibrationPointsDistanceFunction::f(vnl_vector<double> const &x, vnl_vector<double> &fx)
{
  double t_3x, t_3y, t_3z, sz, cz, sy, cy, sx, cx, m_x, m_y;
  double R3_11, R3_21, R3_31, R3_12, R3_22, R3_32;

  t_3x = x(0);
  t_3y = x(1); 
  t_3z = x(2);
  sz = sin(x(3));
  cz = cos(x(3));
  sy = sin(x(4));
  cy = cos(x(4));
  sx = sin(x(5));
  cx = cos(x(5));
  m_x = x(6);
  m_y = x(7);

  R3_11 = cz*cy;
  R3_21 = sz*cy;
  R3_31 = -sy;
  R3_12 = cz*sy*sx - sz*cx;
  R3_22 = sz*sy*sx+cz*cx;
  R3_32 = cy*sx;
 
  unsigned int numPoints = static_cast<unsigned int>(this->data->size());
  double R2[3][3], t2[3];
  double A_11, A_12, A_13, A_14, A_15, A_16, A_17, A_18, A_19;
  double A_21, A_22, A_23, A_24, A_25, A_26, A_27, A_28, A_29;
  double A_31, A_32, A_33, A_34, A_35, A_36, A_37, A_38, A_39;
  double b_1, b_2, b_3, expr1, expr2, expr3;

  for(unsigned int i =0; i<numPoints; i++) { 
    lsqrRecipes::Frame &T2 = (*(this->data))[i]->T2;
    T2.getRotationMatrix(R2);
    T2.getTranslation(t2);
    lsqrRecipes::Point2D &q = (*(this->data))[i]->q;
    lsqrRecipes::Point3D &p = (*(this->data))[i]->p;
             //u_i*R_2i
    A_11 = q[0]*R2[0][0];   A_12 = q[0]*R2[0][1];  A_13 = q[0]*R2[0][2];
    A_21 = q[0]*R2[1][0];   A_22 = q[0]*R2[1][1];  A_23 = q[0]*R2[1][2];
    A_31 = q[0]*R2[2][0];   A_32 = q[0]*R2[2][1];  A_33 = q[0]*R2[2][2];
             //v_i*R_2i
    A_14 = q[1]*R2[0][0];   A_15 = q[1]*R2[0][1];  A_16 = q[1]*R2[0][2];
    A_24 = q[1]*R2[1][0];   A_25 = q[1]*R2[1][1];  A_26 = q[1]*R2[1][2];
    A_34 = q[1]*R2[2][0];   A_35 = q[1]*R2[2][1];  A_36 = q[1]*R2[2][2];
               //R_2i
    A_17 = R2[0][0];   A_18 = R2[0][1];  A_19 = R2[0][2];
    A_27 = R2[1][0];   A_28 = R2[1][1];  A_29 = R2[1][2];
    A_37 = R2[2][0];   A_38 = R2[2][1];  A_39 = R2[2][2];

    b_1 = p[0]-t2[0];
    b_2 = p[1]-t2[1];
    b_3 = p[2]-t2[2];

    expr1 = A_11*m_x*R3_11 +
            A_12*m_x*R3_21 +
            A_13*m_x*R3_31 + 
            A_14*m_y*R3_12 + 
            A_15*m_y*R3_22 + 
            A_16*m_y*R3_32 + 
            A_17*t_3x +
            A_18*t_3y +
            A_19*t_3z -
            b_1;
        
    expr2 = A_21*m_x*R3_11 + 
            A_22*m_x*R3_21 + 
            A_23*m_x*R3_31 + 
            A_24*m_y*R3_12 +
            A_25*m_y*R3_22 +
            A_26*m_y*R3_32 +
            A_27*t_3x +
            A_28*t_3y +
            A_29*t_3z -
            b_2;

    expr3 = A_31*m_x*R3_11 +
            A_32*m_x*R3_21 +
            A_33*m_x*R3_31 +
            A_34*m_y*R3_12 +
            A_35*m_y*R3_22 +
            A_36*m_y*R3_32 +
            A_37*t_3x +
            A_38*t_3y +
            A_39*t_3z -
            b_3;
          //fx has the correct size (see vnl_least_squares_function documentation)
    fx[i] = sqrt(expr1*expr1+expr2*expr2+expr3*expr3);
  }
}

/*****************************************************************************/
void CalibratedPointerTargetUSCalibrationParametersEstimator::SumSquaresCalibrationPointsDistanceFunction::gradf(vnl_vector<double> const& x, vnl_matrix<double>& jacobian)
{
  double t_3x, t_3y, t_3z, sz, cz, sy, cy, sx, cx, m_x, m_y;
  double R3_11, R3_21, R3_31, R3_12, R3_22, R3_32;

  t_3x = x(0);
  t_3y = x(1); 
  t_3z = x(2);
  sz = sin(x(3));
  cz = cos(x(3));
  sy = sin(x(4));
  cy = cos(x(4));
  sx = sin(x(5));
  cx = cos(x(5));
  m_x = x(6);
  m_y = x(7);

  R3_11 = cz*cy;
  R3_21 = sz*cy;
  R3_31 = -sy;
  R3_12 = cz*sy*sx - sz*cx;
  R3_22 = sz*sy*sx+cz*cx;
  R3_32 = cy*sx;

  unsigned int numPoints = static_cast<unsigned int>(this->data->size());
  double R2[3][3], t2[3];
  double A_11, A_12, A_13, A_14, A_15, A_16, A_17, A_18, A_19;
  double A_21, A_22, A_23, A_24, A_25, A_26, A_27, A_28, A_29;
  double A_31, A_32, A_33, A_34, A_35, A_36, A_37, A_38, A_39;
  double b_1, b_2, b_3, expr1, expr2, expr3, delta_i;
  double v1, v2, v3, v4, v5, v6;

  for(unsigned int i =0; i<numPoints; i++) { 
    lsqrRecipes::Frame &T2 = (*(this->data))[i]->T2;
    T2.getRotationMatrix(R2);
    T2.getTranslation(t2);
    lsqrRecipes::Point2D &q = (*(this->data))[i]->q;
    lsqrRecipes::Point3D &p = (*(this->data))[i]->p;
             //u_i*R_2i
    A_11 = q[0]*R2[0][0];   A_12 = q[0]*R2[0][1];  A_13 = q[0]*R2[0][2];
    A_21 = q[0]*R2[1][0];   A_22 = q[0]*R2[1][1];  A_23 = q[0]*R2[1][2];
    A_31 = q[0]*R2[2][0];   A_32 = q[0]*R2[2][1];  A_33 = q[0]*R2[2][2];
             //v_i*R_2i
    A_14 = q[1]*R2[0][0];   A_15 = q[1]*R2[0][1];  A_16 = q[1]*R2[0][2];
    A_24 = q[1]*R2[1][0];   A_25 = q[1]*R2[1][1];  A_26 = q[1]*R2[1][2];
    A_34 = q[1]*R2[2][0];   A_35 = q[1]*R2[2][1];  A_36 = q[1]*R2[2][2];
               //R_2i
    A_17 = R2[0][0];   A_18 = R2[0][1];  A_19 = R2[0][2];
    A_27 = R2[1][0];   A_28 = R2[1][1];  A_29 = R2[1][2];
    A_37 = R2[2][0];   A_38 = R2[2][1];  A_39 = R2[2][2];

    b_1 = p[0]-t2[0];
    b_2 = p[1]-t2[1];
    b_3 = p[2]-t2[2];

    expr1 = A_11*m_x*R3_11 +
            A_12*m_x*R3_21 +
            A_13*m_x*R3_31 + 
            A_14*m_y*R3_12 + 
            A_15*m_y*R3_22 + 
            A_16*m_y*R3_32 + 
            A_17*t_3x +
            A_18*t_3y +
            A_19*t_3z -
            b_1;
        
    expr2 = A_21*m_x*R3_11 + 
            A_22*m_x*R3_21 + 
            A_23*m_x*R3_31 + 
            A_24*m_y*R3_12 +
            A_25*m_y*R3_22 +
            A_26*m_y*R3_32 +
            A_27*t_3x +
            A_28*t_3y +
            A_29*t_3z -
            b_2;

    expr3 = A_31*m_x*R3_11 +
            A_32*m_x*R3_21 +
            A_33*m_x*R3_31 +
            A_34*m_y*R3_12 +
            A_35*m_y*R3_22 +
            A_36*m_y*R3_32 +
            A_37*t_3x +
            A_38*t_3y +
            A_39*t_3z -
            b_3;

    delta_i = sqrt(expr1*expr1+expr2*expr2+expr3*expr3);     
   
             //partial derivative t_3x   
    jacobian(i,0) = (A_17*expr1+A_27*expr2+A_37*expr3)/delta_i;
                //partial derivative t_3y   
    jacobian(i,1) = (A_18*expr1+A_28*expr2+A_38*expr3)/delta_i;
                //partial derivative t_3z   
    jacobian(i,2) = (A_19*expr1+A_29*expr2+A_39*expr3)/delta_i;   
         //partial derivative omega_z
    v1 = -m_x*sz*cy;
    v2 = m_x*cz*cy;
    v3 = -m_y*(sz*sy*sx+cz*cx);
    v4 = m_y*(cz*sy*sx-sz*cx);   
    jacobian(i,3) = ((A_11*v1+A_12*v2+A_14*v3+A_15*v4)*expr1 +
                     (A_21*v1+A_22*v2+A_24*v3+A_25*v4)*expr2 +
                     (A_31*v1+A_32*v2+A_34*v3+A_35*v4)*expr3)/
                     delta_i;
          //partial derivative omega_y
    v1 = -m_x*sy*cz;
    v2 = -m_x*sy*sz;
    v3 = -m_x*cy;
    v4 = m_y*sx*cy*cz;
    v5 = m_y*sx*cy*sz;
    v6 = -m_y*sx*sy;
    jacobian(i,4) = ((A_11*v1+A_12*v2+A_13*v3+A_14*v4+A_15*v5+A_16*v6)*expr1 +
                     (A_21*v1+A_22*v2+A_23*v3+A_24*v4+A_25*v5+A_26*v6)*expr2 +
                     (A_31*v1+A_32*v2+A_33*v3+A_34*v4+A_35*v5+A_36*v6)*expr3)/
                     delta_i;                  
          //partial derivative omega_x
    v1 = m_y*(cz*sy*cx+sz*sx);
    v2 = m_y*(sz*sy*cx-cz*sx);
    v3 = m_y*cy*cx;
    jacobian(i,5) = ((A_14*v1+A_15*v2+A_16*v3)*expr1 +
                     (A_24*v1+A_25*v2+A_26*v3)*expr2 +
                     (A_34*v1+A_35*v2+A_36*v3)*expr3)/
                     delta_i;                  
             //partial derivative m_x
    jacobian(i,6) = ((A_11*R3_11+A_12*R3_21+A_13*R3_31)*expr1 +
                      (A_21*R3_11+A_22*R3_21+A_23*R3_31)*expr2 +
                      (A_31*R3_11+A_32*R3_21+A_33*R3_31)*expr3)/
                      delta_i;
             //partial derivative m_y                   
    jacobian(i,7) = ((A_14*R3_12+A_15*R3_22+A_16*R3_32)*expr1 +
                      (A_24*R3_12+A_25*R3_22+A_26*R3_32)*expr2 +
                      (A_34*R3_12+A_35*R3_22+A_36*R3_32)*expr3)/
                      delta_i;                                      
  }
}


} //namespace lsqrRecipes
