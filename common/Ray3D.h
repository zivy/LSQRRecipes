#ifndef _RAY3D_H_
#define _RAY3D_H_

#include "copyright.h"

#include "Point3D.h"
#include "Vector3D.h"
#include "Frame.h"

/**
 * Class representing 3D rays using a point, p, and a normal, n, such that the ray is r(t) = p+t*n and t is in [0,inf).
 * Note that it is not assumed that ||n|| = 1.
 *
 * @author: Ziv Yaniv (zivy@isis.georgetown.edu)
 */

namespace lsqrRecipes {

class Ray3D {
public:
	friend std::ostream &operator<<(std::ostream& output, const Ray3D &r) { output<<"p: "<<r.p<<" n: "<<r.n; return output;}
	Point3D p;
	Vector3D n;
	Ray3D() {p[0]=p[1]=p[2]=0.0; n[0]=n[1]=n[2]=0;}
  Ray3D(const Ray3D &other) {
    this->n[0] = other.n[0]; this->n[1] = other.n[1]; this->n[2] = other.n[2];
    this->p[0] = other.p[0]; this->p[1] = other.p[1]; this->p[2] = other.p[2];
  }

	void transform(Frame &transformation) {
		transformation.apply(p);
		transformation.apply(n);
	}

 /**
	* This function computes the distance between a point and a ray (line).
	* The Ray is treated as a line , this means that I don't care that
	* the parameter t isn't in [0,infinity).
	*   The parametric equation of the line defined by the ray (Point+direction)
	*   is: 
	*        (*)  P(t) = P0+tN
	*   The closest point on the line P(t) to a given point P3 is on the line that
	*   goes through P3 and is perpendicular to P. 
	*   This gives us: 
	*                (**)  (P3-P(t_i)) dot N = 0;
	*   Substituting (*) into (**) gives:
	*                (x3-x0)*Nx + (y3-y0)*Ny + (z3-z0)*Nz   
	*           t_i= --------------------------------------
	*                        Nx^2 + Ny^2 + Nz^2
	*  
	*   Now we have P(t_i) =P0+ t_iN and the distance between P3 and the
	*   line is the distance between P3 and P(t_i).
	*  
	*/
	double distance(Point3D &pnt);

	/**
   * Compute the intersection of two rays analytically. Based on the code
   * found in Graphics Gems p.304 "Intersection of two 3D lines", Ronald Goldman.
   * I add the definition that the intersection point of two parallel rays is
   * the point between them as shown in the following figure:
   *
   *                 r1(t)*-------------------------->
   *       the intersection-->       X
   *                            r2(s)*--------------->
   * The rays must have a common support. 
   * One of the following is true:
   *   1.Projecting r1(0) onto r2 gives a point with a positive ray parameter value.
   *   2.Projecting r2(0) onto r1 gives a point with a positive ray parameter value.
   *        
  */ 
	bool intersection(Ray3D &r, Point3D &result);
};


/**
 * Lightweight container for ray bundles.
 */
class RayBundle {
public:
	RayBundle(int size) : n(size) {p[0]=p[1]=p[2]=0.0;}
	         //the ray directions (not necessarily normalized)
	std::vector<Vector3D> n;
	        //the rays origin
	Point3D p;

	void transform(Frame &transformation) {
		transformation.apply(p);
		int numRays = static_cast<int>(n.size());
		for(int i=0; i<numRays; i++)
			transformation.apply(n[i]);
	}

	/**
	 * Write the ray bundle to the stream in open inventor ascii format.
	 * The ray bundle is described by a set of lines which originate at the rays origin
	 * and have length 'rayLength'. The origin is marked by a sphere. Both the sphere
	 * and lines are drawn using 'color'.
	 */
	void writeOIVData(std::ostream &out, double rayLength, double color[3]);

};

} //namespace lsqrRecipes

#endif //_RAY3D_H_
