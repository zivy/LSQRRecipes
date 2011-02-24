#include "Ray3D.h"
#include "Epsilon.h"

namespace lsqrRecipes {

bool Ray3D::intersection(Ray3D &r, Point3D &result)
{
	const double INTERSECTION_EPSILON = EPS;
	Vector3D p21, n1Crossn2, tmp, tmp1, tmp2, n1, n2;
	double t, t1, t2;
	
	n1 = this->n;
	n2 = r.n;
	
	p21 = r.p - this->p;

	n1Crossn2 = crossProduct(n1,n2);
	double denominator = n1Crossn2*n1Crossn2;
	
	                 //rays are parallel
	if(denominator < INTERSECTION_EPSILON) {
		if((t=p21*n1) > 0) {
			tmp =  n1*t;
			result[0] = (r.p[0]+ this->p[0]+tmp[0])/2.0;
			result[1] = (r.p[1]+ this->p[1]+tmp[1])/2.0;
			result[2] = (r.p[2]+ this->p[2]+tmp[2])/2.0;
			return true;
		}
		else if((t=-(p21*n2)) > 0) {
			tmp = n2*t;
			result[0] = (r.p[0]+ this->p[0]+tmp[0])/2.0;
			result[1] = (r.p[1]+ this->p[1]+tmp[1])/2.0;
			result[2] = (r.p[2]+ this->p[2]+tmp[2])/2.0;
			return true;
		}
		else return false; //rays are parallel but do not have a common support
	}
      //rays are not parallel (skew or intersecting)
	t1 = (n1Crossn2[0]*(p21[1]*n2[2] - p21[2]*n2[1]) -
			  n1Crossn2[1]*(p21[0]*n2[2] - p21[2]*n2[0]) +
			  n1Crossn2[2]*(p21[0]*n2[1] - p21[1]*n2[0]))/denominator;
	t2 = (n1Crossn2[0]*(p21[1]*n1[2] - p21[2]*n1[1]) -
			  n1Crossn2[1]*(p21[0]*n1[2] - p21[2]*n1[0]) +
			  n1Crossn2[2]*(p21[0]*n1[1] - p21[1]*n1[0]))/denominator;

	if(t1<0 || t2<0)
		return false;

	tmp1 = n1*t1;
	tmp2 = n2*t2;
	
	result[0] = (this->p[0] + tmp1[0] + r.p[0] + tmp2[0])/2.0;
	result[1] = (this->p[1] + tmp1[1] + r.p[1] + tmp2[1])/2.0;
	result[2] = (this->p[2] + tmp1[2] + r.p[2] + tmp2[2])/2.0;
	return true;
}
/*****************************************************************************/
double Ray3D::distance(Point3D &pnt) 
{
  double x0 = this->p[0];
  double y0 = this->p[1];
  double z0 = this->p[2];
  double x3 = pnt[0];
  double y3 = pnt[1];
  double z3 = pnt[2];
  
  double t = this->n[0]*(x3-x0) + this->n[1]*(y3-y0) + this->n[2]*(z3-z0);
  t/= ((this->n[0]*this->n[0]) + (this->n[1]*this->n[1]) + 
       (this->n[2]*this->n[2]));

  double x = x0+t*this->n[0];
  double y = y0+t*this->n[1];
  double z = z0+t*this->n[2];

  return sqrt((x-x3)*(x-x3) + (y-y3)*(y-y3) + (z-z3)*(z-z3));
}
/*****************************************************************************/
void RayBundle::writeOIVData(std::ostream &out, double rayLength, double color[3])
{
	out<<"Separator {\n";
	out<<"\tMaterial {\n";
	out<<"\t\tdiffuseColor "<<color[0]<<" "<<color[1]<<" "<<color[2]<<"\n";
	out<<"\t}\n";
	out<<"\tSeparator {\n";
	out<<"\t\tTransform {\n";
  out<<"\t\t\t translation "<<p[0]<<" "<<p[1]<<" "<<p[2]<<"\n";
	out<<"\t\t}\n";
	out<<"\t\tSphere {\n";
  out<<"\t\t\tradius 10\n";
	out<<"\t\t}\n";
	out<<"\t}\n";
	out<<"\tCoordinate3 {\n";
	out<<"\t\tpoint[\n";
	out<<"\t\t\t"<<this->p[0]<<"\t"<<this->p[1]<<"\t"<<this->p[2]<<",\n";
	int rayNum = static_cast<int>(this->n.size());
	int i;
	for(i=0; i<rayNum; i++) {
		Vector3D &v = this->n[i];
		out<<"\t\t\t"<<this->p[0] + rayLength*v[0]<<"\t"<<this->p[1] + rayLength*v[1]<<"\t"<<this->p[2] + rayLength*v[2]<<",\n";
	}
	out<<"\t\t]\t\n}\n";
	out<<"\tIndexedLineSet {\n";
  out<<"\t\tcoordIndex [\n";
  for(i=0; i<rayNum; i++)
		out<<"\t0, "<<i+1<<", -1,\n";
	out<<"\t]\n}\n}";
}

} //namespace lsqrRecipes
