#include "Frame.h"
#include "AbsoluteOrientationParametersEstimator.h"
#include "Epsilon.h"
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/vnl_cross.h>

namespace lsqrRecipes {

AbsoluteOrientationParametersEstimator::AbsoluteOrientationParametersEstimator(double delta) : ParametersEstimator< std::pair<Point3D,Point3D>,double >(3) {this->deltaSquared = delta*delta;}
/*****************************************************************************/
/*
 * Compute the transformation parameters  [s,q_x,q_y,q_z,t_x,t_y,t_z]
 */
void AbsoluteOrientationParametersEstimator::estimate(std::vector< std::pair<Point3D,Point3D> *> &data, 
                																      std::vector<double> &parameters)
{
	parameters.clear();
	              //not enough data elements for computation
	if(data.size()<this->minForEstimate)
		return;
	vnl_matrix<double> firstR(3,3), secondR(3,3), R(3,3);
  vnl_vector<double> x(3), y(3), z(3), t(3), meanFirst(3), meanSecond(3);

	Point3D &firstP0 = data[0]->first;
	Point3D &firstP1 = data[1]->first;
	Point3D &firstP2 = data[2]->first;
  meanFirst[0] = (firstP0[0] + firstP1[0] + firstP2[0])/3.0;
  meanFirst[1] = (firstP0[1] + firstP1[1] + firstP2[1])/3.0;
  meanFirst[2] = (firstP0[2] + firstP1[2] + firstP2[2])/3.0;

                //define X axis for 'first' set of points
	x[0] = firstP0[0] - meanFirst[0];
	x[1] = firstP0[1] - meanFirst[1];
	x[2] = firstP0[2] - meanFirst[2];
	x.normalize();
               //define Y axis for 'first' set of points
	y[0] = firstP1[0] - meanFirst[0];
	y[1] = firstP1[1] - meanFirst[1];
	y[2] = firstP1[2] - meanFirst[2];
	y = y - dot_product(y,x)*x;
	y.normalize();
              //define Z axis for 'first' set of points
	z = vnl_cross_3d(x,y);
                       
  if(z.magnitude()<EPS) //points are collinear
    return;

	                //place axes in rotation matrix
	firstR.set_column(0,x);
	firstR.set_column(1,y);
	firstR.set_column(2,z);

	Point3D &secondP0 = data[0]->second;
	Point3D &secondP1 = data[1]->second;
	Point3D &secondP2 = data[2]->second;
  meanSecond[0] = (secondP0[0] + secondP1[0] + secondP2[0])/3.0;
  meanSecond[1] = (secondP0[1] + secondP1[1] + secondP2[1])/3.0;
  meanSecond[2] = (secondP0[2] + secondP1[2] + secondP2[2])/3.0;

                //define X axis for 'second' set of points
	x[0] = secondP0[0] - meanSecond[0];
	x[1] = secondP0[1] - meanSecond[1];
	x[2] = secondP0[2] - meanSecond[2];
	x.normalize();
               //define Y axis for 'second' set of points
	y[0] = secondP1[0] - meanSecond[0];
	y[1] = secondP1[1] - meanSecond[1];
	y[2] = secondP1[2] - meanSecond[2];
	y = y - dot_product(y,x)*x;
	y.normalize();
              //define Z axis for 'second' set of points
	z = vnl_cross_3d(x,y);

                        //if we got here then this shouldn't happen (better safe than sorry)
  if(z.magnitude()<EPS) //points are collinear
    return;

          //place axes in rotation matrix
	secondR.set_column(0,x);
	secondR.set_column(1,y);
	secondR.set_column(2,z);

	                //compute rotation
	R = secondR*firstR.transpose();
                 //compute translation   
  t = meanSecond - R*meanFirst;

	Frame frm;
	frm.setRotationMatrix(R(0,0), R(0,1), R(0,2),
												R(1,0), R(1,1), R(1,2),
                        R(2,0), R(2,1), R(2,2));
	double quaternion[4];
	frm.getRotationQuaternion(quaternion);
  parameters.push_back(quaternion[0]);
	parameters.push_back(quaternion[1]);
	parameters.push_back(quaternion[2]);
	parameters.push_back(quaternion[3]);
	parameters.push_back(t[0]);
	parameters.push_back(t[1]);
	parameters.push_back(t[2]);
}
/*****************************************************************************/
/*
 * Compute the transformation parameters  [s,q_x,q_y,q_z,t_x,t_y,t_z]
 */
void AbsoluteOrientationParametersEstimator::estimate(std::vector< std::pair<Point3D,Point3D> > &data, 
																	                    std::vector<double> &parameters)
{
	std::vector< std::pair<Point3D,Point3D> *> usedData;
	int dataSize = data.size();
	for(int i=0; i<dataSize; i++)
		usedData.push_back(&(data[i]));
	estimate(usedData,parameters);
}
/*****************************************************************************/
/*
 * Compute the transformation parameters  [s,q_x,q_y,q_z,t_x,t_y,t_z]
 * NEED TO ADD CHECKING FOR DEGENERATE CONFIGURATION.
 */
void AbsoluteOrientationParametersEstimator::leastSquaresEstimate(std::vector< std::pair<Point3D,Point3D> *> &data, 
																							                    std::vector<double> &parameters)
{
	parameters.clear();
	unsigned int i, pairNum = data.size();
	              //not enough data elements for computation
	if(pairNum<this->minForEstimate)
		return;

	vnl_matrix<double> muLmuR(3,3,0), M(3,3,0), curMat(3,3,0), N(4,4,0);
	Point3D meanFirst, meanSecond; //assume points set to zero by constructor
	

	       //compute the mean of both point sets
	for (i=0; i<pairNum; i++) {
		meanFirst[0] += data[i]->first[0];	    meanFirst[1] += data[i]->first[1];	    meanFirst[2] += data[i]->first[2];
		meanSecond[0] += data[i]->second[0];	  meanSecond[1] += data[i]->second[1];	  meanSecond[2] += data[i]->second[2];
	}
	meanFirst[0]/=pairNum;	  meanFirst[1]/=pairNum;	  meanFirst[2]/=pairNum;
	meanSecond[0]/=pairNum;	  meanSecond[1]/=pairNum;	  meanSecond[2]/=pairNum;

              //compute the matrix muLmuR
	muLmuR(0,0) = meanFirst[0]*meanSecond[0];		
	muLmuR(0,1) = meanFirst[0]*meanSecond[1];		
	muLmuR(0,2) = meanFirst[0]*meanSecond[2];
	muLmuR(1,0) = meanFirst[1]*meanSecond[0];
	muLmuR(1,1) = meanFirst[1]*meanSecond[1];
	muLmuR(1,2) = meanFirst[1]*meanSecond[2];
	muLmuR(2,0) = meanFirst[2]*meanSecond[0];
	muLmuR(2,1) = meanFirst[2]*meanSecond[1];
	muLmuR(2,2) = meanFirst[2]*meanSecond[2];

	                    //compute the matrix M
	for (i=0; i<pairNum; i++) {
		Point3D &leftPoint = data[i]->first;
		Point3D &rightPoint = data[i]->second;
		curMat(0,0) = leftPoint[0]*rightPoint[0];		
		curMat(0,1) = leftPoint[0]*rightPoint[1];		
		curMat(0,2) = leftPoint[0]*rightPoint[2];
		curMat(1,0) = leftPoint[1]*rightPoint[0];
		curMat(1,1) = leftPoint[1]*rightPoint[1];
		curMat(1,2) = leftPoint[1]*rightPoint[2];
		curMat(2,0) = leftPoint[2]*rightPoint[0];
		curMat(2,1) = leftPoint[2]*rightPoint[1];
		curMat(2,2) = leftPoint[2]*rightPoint[2];
		M+=curMat;
	}

	M+= (muLmuR* -(static_cast<int>(pairNum)));

            	//compute the matrix N	
	vnl_matrix<double> tmpMat(3,3,0);
	double A12, A20, A01;
  double traceM = 0.0;
  for(i=0; i<3; i++)
    traceM+=M(i,i);

	tmpMat.fill_diagonal(-traceM);
	tmpMat += (M + M.transpose());

  A12 = M(1,2) - M(2,1);
  A20 = M(2,0) - M(0,2);
  A01 = M(0,1) - M(1,0);

  N(0,0)=traceM; N(0,1)=A12; N(0,2)=A20; N(0,3)=A01;
  N(1,0)=A12;
  N(2,0)=A20;
  N(3,0)=A01;
  N.update(tmpMat,1,1);

            //find the eigenvector that belongs to the maximal 
           //eigenvalue of N, eigenvalues are sorted from smallest to largest
	vnl_symmetric_eigensystem<double> eigenSystem(N);
                          
	   //the transformation parameters  [s,q_x,q_y,q_z,t_x,t_y,t_z]
	parameters.push_back(eigenSystem.V(0,3));
	parameters.push_back(eigenSystem.V(1,3));
	parameters.push_back(eigenSystem.V(2,3));
	parameters.push_back(eigenSystem.V(3,3));

	Frame frm;
	frm.setRotationQuaternion(eigenSystem.V(0,3),eigenSystem.V(1,3),eigenSystem.V(2,3),eigenSystem.V(3,3), true);
	frm.apply(meanFirst);	
	parameters.push_back(meanSecond[0] - meanFirst[0]);
	parameters.push_back(meanSecond[1] - meanFirst[1]);
	parameters.push_back(meanSecond[2] - meanFirst[2]);
}
/*****************************************************************************/
void AbsoluteOrientationParametersEstimator::weightedLeastSquaresEstimate(std::vector< std::pair<Point3D,Point3D> *> &data, std::vector<double> &weights,
		                                                                      std::vector<double> &parameters)
{
	parameters.clear();
	unsigned int i, pairNum = data.size();
	              //not enough data elements for computation
	if(pairNum<this->minForEstimate)
		return;

	vnl_matrix<double> muLmuR(3,3,0), M(3,3,0), curMat(3,3,0), N(4,4,0);
	Point3D meanFirst, meanSecond; //assume points set to zero by constructor
	double sumWeights;

	sumWeights = 0.0;
	for(i=0; i<pairNum; i++)
		sumWeights += weights[i];

	       //compute the mean of both point sets
	for (i=0; i<pairNum; i++) {
		meanFirst[0] += data[i]->first[0]*weights[i];	    meanFirst[1] += data[i]->first[1]*weights[i];	    meanFirst[2] += data[i]->first[2]*weights[i];
		meanSecond[0] += data[i]->second[0]*weights[i];	  meanSecond[1] += data[i]->second[1]*weights[i];	  meanSecond[2] += data[i]->second[2]*weights[i];
	}
	meanFirst[0]/=sumWeights;	    meanFirst[1]/=sumWeights;	    meanFirst[2]/=sumWeights;
	meanSecond[0]/=sumWeights;	  meanSecond[1]/=sumWeights;	  meanSecond[2]/=sumWeights;

	              //compute the matrix muLmuR
	muLmuR(0,0) = meanFirst[0]*meanSecond[0];		
	muLmuR(0,1) = meanFirst[0]*meanSecond[1];		
	muLmuR(0,2) = meanFirst[0]*meanSecond[2];
	muLmuR(1,0) = meanFirst[1]*meanSecond[0];
	muLmuR(1,1) = meanFirst[1]*meanSecond[1];
	muLmuR(1,2) = meanFirst[1]*meanSecond[2];
	muLmuR(2,0) = meanFirst[2]*meanSecond[0];
	muLmuR(2,1) = meanFirst[2]*meanSecond[1];
	muLmuR(2,2) = meanFirst[2]*meanSecond[2];

	                    //compute the matrix M
	for (i=0; i<pairNum; i++) {
		Point3D &leftPoint = data[i]->first;
		Point3D &rightPoint = data[i]->second;
		curMat(0,0) = leftPoint[0]*rightPoint[0];		
		curMat(0,1) = leftPoint[0]*rightPoint[1];		
		curMat(0,2) = leftPoint[0]*rightPoint[2];
		curMat(1,0) = leftPoint[1]*rightPoint[0];
		curMat(1,1) = leftPoint[1]*rightPoint[1];
		curMat(1,2) = leftPoint[1]*rightPoint[2];
		curMat(2,0) = leftPoint[2]*rightPoint[0];
		curMat(2,1) = leftPoint[2]*rightPoint[1];
		curMat(2,2) = leftPoint[2]*rightPoint[2];
		M+=(curMat*weights[i]);
	}
	M+= (muLmuR *(-sumWeights));

            	//compute the matrix N	
	vnl_matrix<double> tmpMat(3,3,0);
	double A12, A20, A01;
  double traceM = 0.0;
  for(i=0; i<3; i++)
    traceM+=M(i,i);

	tmpMat.fill_diagonal(-traceM);
	tmpMat += (M + M.transpose());

  A12 = M(1,2) - M(2,1);
  A20 = M(2,0) - M(0,2);
  A01 = M(0,1) - M(1,0);

  N(0,0)=traceM; N(0,1)=A12; N(0,2)=A20; N(0,3)=A01;
  N(1,0)=A12;
  N(2,0)=A20;
  N(3,0)=A01;
  N.update(tmpMat,1,1);

            //find the eigenvector that belongs to the maximal 
           //eigenvalue of N, eigenvalues are sorted from smallest to largest
	vnl_symmetric_eigensystem<double> eigenSystem(N);
                          
	   //the transformation parameters  [s,q_x,q_y,q_z,t_x,t_y,t_z]
	parameters.push_back(eigenSystem.V(0,3));
	parameters.push_back(eigenSystem.V(1,3));
	parameters.push_back(eigenSystem.V(2,3));
	parameters.push_back(eigenSystem.V(3,3));

	Frame frm;
	frm.setRotationQuaternion(eigenSystem.V(0,3),eigenSystem.V(1,3),eigenSystem.V(2,3),eigenSystem.V(3,3), true);
	frm.apply(meanFirst);	
	parameters.push_back(meanSecond[0] - meanFirst[0]);
	parameters.push_back(meanSecond[1] - meanFirst[1]);
	parameters.push_back(meanSecond[2] - meanFirst[2]);
}
/*****************************************************************************/
/*
 * Compute the transformation parameters  [s,q_x,q_y,q_z,t_x,t_y,t_z]
 */
void AbsoluteOrientationParametersEstimator::leastSquaresEstimate(std::vector< std::pair<Point3D,Point3D> > &data, 
																							                    std::vector<double> &parameters)
{
	std::vector< std::pair<Point3D,Point3D> *> usedData;
	int dataSize = data.size();
	for(int i=0; i<dataSize; i++)
		usedData.push_back(&(data[i]));
	leastSquaresEstimate(usedData,parameters);
}
/*****************************************************************************/
/*
 * Given the transformation parameters  [s,q_x,q_y,q_z,t_x,t_y,t_z] check if
 * ||(data.second - T(parameters)*data.first)||^2 < delta^2
 */
bool AbsoluteOrientationParametersEstimator::agree(std::vector<double> &parameters, std::pair<Point3D,Point3D> &data)
{
	Frame transformation(parameters[4], parameters[5], parameters[6], 
		                   parameters[0], parameters[1], parameters[2], parameters[3]);
	Point3D p1Transformed;
	transformation.apply(data.first, p1Transformed);

	double dx = p1Transformed[0] - data.second[0];
	double dy = p1Transformed[1] - data.second[1];
	double dz = p1Transformed[2] - data.second[2];
	return ((dx*dx + dy*dy + dz*dz) < this->deltaSquared);
}

} //namespace lsqrRecipes
