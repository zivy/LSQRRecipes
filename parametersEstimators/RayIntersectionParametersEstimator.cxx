#include "RayIntersectionParametersEstimator.h"
#include <vnl/algo/vnl_matrix_inverse.h>
#include "Epsilon.h"
#include "Vector3D.h"


namespace lsqrRecipes {

RayIntersectionParametersEstimator::RayIntersectionParametersEstimator(double delta, 
                                                                       double minimalAngularDeviation) : 
  ParametersEstimator<Ray3D,double>(2) 
{
  this->deltaSquared = delta*delta; 
  this->crossEps = sin(minimalAngularDeviation); 
  this->crossEps*=this->crossEps;
}
/*****************************************************************************/
/**
 * Compute the intersection of two rays. Based on the code
 * found in Graphics Gems p.304 "Intersection of Two Lines in Three-Space", 
 * Ronald Goldman.        
 */ 
void RayIntersectionParametersEstimator::estimate(std::vector<Ray3D *> &data, 
																	                std::vector<double> &parameters)
{
	parameters.clear();
	              //not enough data elements for computation
	if(data.size()<this->minForEstimate)
		return;

  lsqrRecipes::Vector3D &n1 =  data[0]->n;
  lsqrRecipes::Vector3D &n2 = data[1]->n;
  lsqrRecipes::Point3D &p1 = data[0]->p;
  lsqrRecipes::Point3D &p2 = data[1]->p;

  lsqrRecipes::Vector3D p21, n1Crossn2;
  double t1, t2;
  	
  p21[0] = p2[0] - p1[0];
  p21[1] = p2[1] - p1[1];
  p21[2] = p2[2] - p1[2];

  n1Crossn2[0] = n1[1]*n2[2] - n1[2]*n2[1];
  n1Crossn2[1] = n1[2]*n2[0] - n1[0]*n2[2];
  n1Crossn2[2] = n1[0]*n2[1] - n1[1]*n2[0];

  double denominator = n1Crossn2[0]*n1Crossn2[0] + 
                       n1Crossn2[1]*n1Crossn2[1] + 
                       n1Crossn2[2]*n1Crossn2[2];
	
                   //rays are parallel they don't intersect (this is where I 
                   //make the assumption that ||n1||=||n2||=1 as the condition 
                   //is ||cross(n1,n2)||^2<(||n1||^2*||n2||^2sin(minimalAngularDeviation)^2 = sin(minimalAngularDeviation)^2 = this->crossEps
  if(denominator<this->crossEps)
    return;
                 //rays are not parallel (skew or intersecting)
	t1 = (n1Crossn2[0]*(p21[1]*n2[2] - p21[2]*n2[1]) -
			  n1Crossn2[1]*(p21[0]*n2[2] - p21[2]*n2[0]) +
			  n1Crossn2[2]*(p21[0]*n2[1] - p21[1]*n2[0]))/denominator;
	t2 = (n1Crossn2[0]*(p21[1]*n1[2] - p21[2]*n1[1]) -
			  n1Crossn2[1]*(p21[0]*n1[2] - p21[2]*n1[0]) +
			  n1Crossn2[2]*(p21[0]*n1[1] - p21[1]*n1[0]))/denominator;
   //the lines {t in (-inf, inf)} intersect, but the rays {t in [0, inf)} don't
	if(t1<0 || t2<0)
		return;
                         //the intersection is the mid-point  
  parameters.push_back((p1[0]+t1*n1[0] + p2[0]+t2*n2[0])/2.0);
  parameters.push_back((p1[1]+t1*n1[1] + p2[1]+t2*n2[1])/2.0);
  parameters.push_back((p1[2]+t1*n1[2] + p2[2]+t2*n2[2])/2.0);
}
/*****************************************************************************/
/*
* Compute the intersection point parameters  [x,y,z]
 */
void RayIntersectionParametersEstimator::estimate(std::vector<Ray3D> &data, 
																	                std::vector<double> &parameters)
{
	std::vector<Ray3D *> usedData;
	size_t dataSize = data.size();
	for(size_t i=0; i<dataSize; i++)
		usedData.push_back(&(data[i]));
	estimate(usedData,parameters);
}
/*****************************************************************************/
/*
 * Given a set of rays which are assumed to intersect at a common point we
 * compute the least squares estimate of this point.
 * A ray is defined as the linear entity: r(t) = p + t*n , t in [0,infinity) .
 *
 * The least squares solution is obtained by solving the following set of 
 * equations:
 *  m - number of rays.
 *  n_i - ray direction, column vector.
 *  p_i - ray origin (the point where the parameter t=0), column vector.
 * 
 *  [diag(m) - sum(n_i*n_i^T)]x = sum(p_i - n_i^Tp_in_i)
 *             ||                          ||
 *              A             x =           b
 */
void RayIntersectionParametersEstimator::leastSquaresEstimate(std::vector<Ray3D *> &data, 
																							                std::vector<double> &parameters)
{
	vnl_matrix<double> A(3,3,0);
	vnl_vector<double> b(3,0), x(3);
	size_t rayNum = data.size();

	          //create the matrix A and vector b
	for(size_t i=0; i<rayNum; i++) {
		Ray3D &ray = *(data[i]);
		A(0,0) += -(ray.n[0] * ray.n[0]);
		A(0,1) += -(ray.n[0] * ray.n[1]);
		A(0,2) += -(ray.n[0] * ray.n[2]);
		A(1,1) += -(ray.n[1] * ray.n[1]);
		A(1,2) += -(ray.n[1] * ray.n[2]);
		A(2,2) += -(ray.n[2] * ray.n[2]);

		double s = ray.n[0] * ray.p[0] + 
							 ray.n[1] * ray.p[1] + 
							 ray.n[2] * ray.p[2];
		
		b[0] += ray.p[0] - s*ray.n[0];
		b[1] += ray.p[1] - s*ray.n[1];
		b[2] += ray.p[2] - s*ray.n[2];
	}
  A(0,0) += rayNum;
  A(1,0) = A(0,1);
	A(1,1) += rayNum;
  A(2,0) = A(0,2);
  A(2,1) = A(1,2);
	A(2,2) += rayNum;

	vnl_matrix_inverse<double> Ainv(A);
               //explicitly zero out small singular values 
               //this is ugly as it exposes that the inverse is computed via SVD
  Ainv.zero_out_absolute(EPS);

	if(Ainv.rank()<3) //all rays are parallel
		return;
	x = Ainv * b;

  parameters.push_back(x[0]);
  parameters.push_back(x[1]);
  parameters.push_back(x[2]);
}
/*****************************************************************************/
/*
 * Compute the intersection point parameters [x,y,z]
 */
void RayIntersectionParametersEstimator::leastSquaresEstimate(std::vector<Ray3D> &data, 
																							std::vector<double> &parameters)
{
	std::vector<Ray3D *> usedData;
	size_t dataSize = data.size();
	for(size_t i=0; i<dataSize; i++)
		usedData.push_back(&(data[i]));
	leastSquaresEstimate(usedData,parameters);
}
/*****************************************************************************/
/*
 * Given the intersection point [x,y,z] check if it is coincident with the ray.
 * If the parameters vector contains less than three entries an exception will
 * be thrown by the vector's [] operator. We don't check vector length explicitly.
 */
bool RayIntersectionParametersEstimator::agree(std::vector<double> &parameters, 
                                               Ray3D &data)
{
  lsqrRecipes::Vector3D &n = data.n;
  lsqrRecipes::Point3D &p = data.p;
  
  double t = n[0]*(parameters[0]-p[0]) + 
             n[1]*(parameters[1]-p[1]) + 
             n[2]*(parameters[2]-p[2]);
  double dx = parameters[0] - p[0] - t*n[0];
  double dy = parameters[1] - p[1] - t*n[1];
  double dz = parameters[2] - p[2] - t*n[2];

              //closest point is on the ray and the distance is less than delta
  return t>=0 && (dx*dx+dy*dy+dz*dz < this->deltaSquared);
}

} //namespace lsqrRecipes
