#include <memory.h>
#include <math.h>
#include "Frame.h"

namespace lsqrRecipes {

                   //0.5 degrees
const double Frame::SMALL_ANGLE = 0.008726535498373935;
const double Frame::PI = 3.14159265358979323846;
const double Frame::HALF_PI = 3.14159265358979323846/2.0;
const double Frame::DEGREES_PER_RADIAN = 180.0/3.14159265358979323846;
const double Frame::RADIANS_PER_DEGREE = 3.14159265358979323846/180.0;

std::ostream &operator<<(std::ostream& output, const Frame &f)
{
	output<<"translation:\n"<<"\t[ "<<f.translation[0]<<", "<<
		f.translation[1]<<", "<<f.translation[2]<<"]\n";
	switch(f.outputFormat) {
	case Frame::EULER_ANGLES:
		double eulerAngles[6];
		bool isGimbalLock;
		f.getRotationAngles(eulerAngles,isGimbalLock);
		output<<"rotation ([ax,ay,az])";
		if(isGimbalLock)
			output<<" HAS GIMBAL LOCK";
		output<<":\n";
		output<<"\t["<<Frame::toDegrees(eulerAngles[0])<<", "<<
			           Frame::toDegrees(eulerAngles[1])<<", "<<
					   Frame::toDegrees(eulerAngles[2])<<"]\n";
		break;
	case Frame::AXIS_ANGLE:
		double axisAndAngle[4];
		f.getRotationAxisAngle(axisAndAngle);
		output<<"rotation ([axis], angle):\n";
		output<<"\t["<<axisAndAngle[1]<<", "<<axisAndAngle[2]<<", "<<axisAndAngle[3]<<"]\n";
		output<<"\t"<<Frame::toDegrees(axisAndAngle[0])<<"\n";
		break;
	case Frame::QUATERNION:
		double quaternion[4];
		f.getRotationQuaternion(quaternion);
		output<<"rotation ([s,qx,qy,qz]):\n";
		output<<"\t["<<quaternion[0]<<", "<<quaternion[1]<<", "<<quaternion[2]<<", "<<
			quaternion[3]<<"]\n";
		break;
	case Frame::MATRIX:
	default:
		output<<"rotation:\n"<<"\t["<<f.rotation[0][0]<<", "<<f.rotation[0][1]<<", "<<f.rotation[0][2]<<"]\n"<<
			    "\t["<<f.rotation[1][0]<<", "<<f.rotation[1][1]<<", "<<f.rotation[1][2]<<"]\n"<<
				"\t["<<f.rotation[2][0]<<", "<<f.rotation[2][1]<<", "<<f.rotation[2][2]<<"]\n";
		break;
	}
  return(output);
}
/*****************************************************************************/
Frame::Frame(double R[3][3], double t[3]) : outputFormat(Frame::MATRIX)
{
	int sz = 3*sizeof(double);
	if(t == NULL) 
		memset(this->translation,'\0',sz);
	else
		memcpy(this->translation, t, sz);

	if(R == NULL) {
		memset(&(this->rotation[0]),'\0',sz);
		memset(&(this->rotation[1]),'\0',sz);
		memset(&(this->rotation[2]),'\0',sz);
		this->rotation[0][0] = 1;
		this->rotation[1][1] = 1;
		this->rotation[2][2] = 1;
	}
	else {
		memcpy(&(this->rotation[0]),&(R[0]),sz);
		memcpy(&(this->rotation[1]),&(R[1]),sz);
		memcpy(&(this->rotation[2]),&(R[2]),sz);
	}
}
/*****************************************************************************/
Frame::Frame(const Frame &f) : outputFormat(Frame::MATRIX)
{
	int sz = 3*sizeof(double);
	memcpy(this->translation,f.translation,sz);
	memcpy(&(this->rotation[0]),&(f.rotation[0]),sz);
	memcpy(&(this->rotation[1]),&(f.rotation[1]),sz);
	memcpy(&(this->rotation[2]),&(f.rotation[2]),sz);
}
/*****************************************************************************/
Frame::Frame(double x, double y, double z, double ax, double ay, double az) : outputFormat(Frame::MATRIX)
{
	double cx, cy, cz, sx, sy, sz;

	this->translation[0] = x;
	this->translation[1] = y;
	this->translation[2] = z;

	cx = cos(ax);
	cy = cos(ay);
	cz = cos(az);
	sx = sin(ax);
	sy = sin(ay);
	sz = sin(az);

	this->rotation[0][0] = cz*cy; 
	this->rotation[0][1] = cz*sy*sx - sz*cx;  
	this->rotation[0][2] = cz*sy*cx+sz*sx;     

	this->rotation[1][0] = sz*cy; 
	this->rotation[1][1] = sz*sy*sx + cz*cx ; 
	this->rotation[1][2] = sz*sy*cx - cz*sx;   

	this->rotation[2][0] = -sy;   
	this->rotation[2][1] = cy*sx;             
	this->rotation[2][2] = cy*cx;
}
/*****************************************************************************/
Frame::Frame(double x, double y, double z, double axis[3], double angle) : outputFormat(Frame::MATRIX)
{
	double c,s,v;

	this->translation[0] = x;
	this->translation[1] = y;
	this->translation[2] = z;


	c =  cos(angle);
	s =  sin(angle);
	v = 1-c;

	this->rotation[0][0] = axis[0]*axis[0]*v + c; 
	this->rotation[0][1] = axis[0]*axis[1]*v - axis[2]*s;  
	this->rotation[0][2] = axis[0]*axis[2]*v + axis[1]*s;     

	this->rotation[1][0] = axis[0]*axis[1]*v + axis[2]*s; 
	this->rotation[1][1] = axis[1]*axis[1]*v + c; 
	this->rotation[1][2] = axis[1]*axis[2]*v - axis[0]*s;   
	
	this->rotation[2][0] = axis[0]*axis[2]*v - axis[1]*s;   
	this->rotation[2][1] = axis[1]*axis[2]*v + axis[0]*s;             
	this->rotation[2][2] = axis[2]*axis[2]*v + c;
}
/*****************************************************************************/
Frame::Frame(double x, double y, double z, double axisAngle[3]) : outputFormat(Frame::MATRIX)
{
	double angle,c,s,v;
	double axis[3];

	this->translation[0] = x;
	this->translation[1] = y;
	this->translation[2] = z;

	angle = sqrt(axisAngle[0]*axisAngle[0] +
         		 axisAngle[1]*axisAngle[1] +
		         axisAngle[2]*axisAngle[2]);
	axis[0] = axisAngle[0]/angle;
	axis[1] = axisAngle[1]/angle;
	axis[2] = axisAngle[2]/angle;

	c =  cos(angle);
	s =  sin(angle);
	v = 1-c;

	this->rotation[0][0] = axis[0]*axis[0]*v + c; 
	this->rotation[0][1] = axis[0]*axis[1]*v - axis[2]*s;  
	this->rotation[0][2] = axis[0]*axis[2]*v + axis[1]*s;     

	this->rotation[1][0] = axis[0]*axis[1]*v + axis[2]*s; 
	this->rotation[1][1] = axis[1]*axis[1]*v + c; 
	this->rotation[1][2] = axis[1]*axis[2]*v - axis[0]*s;   
	
	this->rotation[2][0] = axis[0]*axis[2]*v - axis[1]*s;   
	this->rotation[2][1] = axis[1]*axis[2]*v + axis[0]*s;             
	this->rotation[2][2] = axis[2]*axis[2]*v + c;
}
/*****************************************************************************/
Frame::Frame(double x, double y, double z, double s, double qx, double qy, double qz, bool normalizeQuaternion) : outputFormat(Frame::MATRIX)
{
	this->translation[0] = x;
	this->translation[1] = y;
	this->translation[2] = z;

	if(normalizeQuaternion) {
		double norm = sqrt(s*s + qx*qx + qy*qy + qz*qz);
		s/=norm;
		qx/=norm;
		qy/=norm;
		qz/=norm;
	}

	this->rotation[0][0] = 1-2*(qy*qy+qz*qz);   
	this->rotation[0][1] = 2*(qx*qy-s*qz);    
	this->rotation[0][2] = 2*(qx*qz+s*qy);    
	
	this->rotation[1][0] = 2*(qx*qy+s*qz);     
	this->rotation[1][1] = 1-2*(qx*qx+qz*qz);  
	this->rotation[1][2] = 2*(qy*qz-s*qx);    
	
	this->rotation[2][0] = 2*(qx*qz-s*qy);     
	this->rotation[2][1] = 2*(qy*qz+s*qx);    
	this->rotation[2][2] = 1-2*(qx*qx+qy*qy);  
}
/*****************************************************************************/
int Frame::setOutputFormat(int newFormat)
{
	int oldFormat = this->outputFormat;
	this->outputFormat = newFormat;
	return oldFormat;
}
/*****************************************************************************/
void Frame::apply(Point3D &p) const
{
	double x, y, z;

	x = this->rotation[0][0] * p[0] + 
	    this->rotation[0][1] * p[1] + 
	    this->rotation[0][2] * p[2] + 
	    this->translation[0];
	y = this->rotation[1][0] * p[0] + 
	    this->rotation[1][1] * p[1] + 
	    this->rotation[1][2] * p[2] + 
	    this->translation[1];
	z = this->rotation[2][0] * p[0] + 
	    this->rotation[2][1] * p[1] + 
	    this->rotation[2][2] * p[2] + 
	    this->translation[2];
	p[0] = x;
	p[1] = y;
	p[2] = z;
}
/*****************************************************************************/
void Frame::apply(const Point3D &p, Point3D &pTransformed) const
{
	double x, y, z;

	x = this->rotation[0][0] * p[0] + 
	    this->rotation[0][1] * p[1] + 
	    this->rotation[0][2] * p[2] + 
	    this->translation[0];
	y = this->rotation[1][0] * p[0] + 
	    this->rotation[1][1] * p[1] + 
	    this->rotation[1][2] * p[2] + 
	    this->translation[1];
	z = this->rotation[2][0] * p[0] + 
	    this->rotation[2][1] * p[1] + 
	    this->rotation[2][2] * p[2] + 
	    this->translation[2];
	pTransformed[0] = x;
	pTransformed[1] = y;
	pTransformed[2] = z;
}
/*****************************************************************************/
void Frame::applyInverse(Point3D &p) 
{
	Frame f;
	f.invert(*this);
	f.apply(p);
}
/*****************************************************************************/
void Frame::applyInverse(Point3D &p, Point3D &pTransformed) 
{
	Frame f;
	f.invert(*this);
	f.apply(p,pTransformed);
}
/*****************************************************************************/
void Frame::apply(Vector3D &v) const
{
	double x, y, z;

	x = this->rotation[0][0] * v[0] + 
	    this->rotation[0][1] * v[1] + 
	    this->rotation[0][2] * v[2];
	y = this->rotation[1][0] * v[0] + 
	    this->rotation[1][1] * v[1] + 
	    this->rotation[1][2] * v[2];
	z = this->rotation[2][0] * v[0] + 
	    this->rotation[2][1] * v[1] + 
	    this->rotation[2][2] * v[2];
	v[0] = x;
	v[1] = y;
	v[2] = z;
}
/*****************************************************************************/
#ifdef USING_VNL
void Frame::apply(vnl_vector_fixed<double,3> &v) const
{
	double x, y, z;

	x = this->rotation[0][0] * v[0] + 
	    this->rotation[0][1] * v[1] + 
	    this->rotation[0][2] * v[2];
	y = this->rotation[1][0] * v[0] + 
	    this->rotation[1][1] * v[1] + 
	    this->rotation[1][2] * v[2];
	z = this->rotation[2][0] * v[0] + 
	    this->rotation[2][1] * v[1] + 
	    this->rotation[2][2] * v[2];
	v[0] = x;
	v[1] = y;
	v[2] = z;
}
#endif
/*****************************************************************************/
void Frame::apply(Vector3D &v, Vector3D &vTransformed) const
{
	double x, y, z;

	x = this->rotation[0][0] * v[0] + 
	    this->rotation[0][1] * v[1] + 
	    this->rotation[0][2] * v[2];
	y = this->rotation[1][0] * v[0] + 
	    this->rotation[1][1] * v[1] + 
	    this->rotation[1][2] * v[2];
	z = this->rotation[2][0] * v[0] + 
	    this->rotation[2][1] * v[1] + 
	    this->rotation[2][2] * v[2];
	vTransformed[0] = x;
	vTransformed[1] = y;
	vTransformed[2] = z;
}
/*****************************************************************************/
#ifdef USING_VNL
void Frame::apply(vnl_vector_fixed<double,3> &v, vnl_vector_fixed<double,3> &vTransformed) const
{
	double x, y, z;

	x = this->rotation[0][0] * v[0] + 
	    this->rotation[0][1] * v[1] + 
	    this->rotation[0][2] * v[2];
	y = this->rotation[1][0] * v[0] + 
	    this->rotation[1][1] * v[1] + 
	    this->rotation[1][2] * v[2];
	z = this->rotation[2][0] * v[0] + 
	    this->rotation[2][1] * v[1] + 
	    this->rotation[2][2] * v[2];
	vTransformed[0] = x;
	vTransformed[1] = y;
	vTransformed[2] = z;
}
#endif
/*****************************************************************************/
void Frame::applyInverse(Vector3D &v) 
{
	Frame f;
	f.invert(*this);
	f.apply(v);
}
/*****************************************************************************/
#ifdef USING_VNL
void Frame::applyInverse(vnl_vector_fixed<double,3> &v) 
{
	Frame f;
	f.invert(*this);
	f.apply(v);
}
#endif
/*****************************************************************************/
void Frame::applyInverse(Vector3D &v, Vector3D &vTransformed) 
{
	Frame f;
	f.invert(*this);
	f.apply(v,vTransformed);
}
/*****************************************************************************/
#ifdef USING_VNL
void Frame::applyInverse(vnl_vector_fixed<double,3> &v, vnl_vector_fixed<double,3> &vTransformed) 
{
	Frame f;
	f.invert(*this);
	f.apply(v,vTransformed);
}
#endif
/*****************************************************************************/
void Frame::mul(const Frame &f1, const Frame &f2) 
{
	int i,j;
  double tmpRotation[3][3];
  double tmpTranslation[3];

  for (i=0;i<3;i++) {
		for (j=0;j<3;j++) {
			tmpRotation[i][j]=f1.rotation[i][0]*f2.rotation[0][j] +
	                      f1.rotation[i][1]*f2.rotation[1][j] +
	                      f1.rotation[i][2]*f2.rotation[2][j];
		}
  }
  for(i=0; i<3; i++) 
		tmpTranslation[i] = f1.rotation[i][0]*f2.translation[0] +
	                      f1.rotation[i][1]*f2.translation[1] +
	                      f1.rotation[i][2]*f2.translation[2] +
	                      f1.translation[i];

	int sz = 3*sizeof(double);
	memcpy(this->translation,tmpTranslation,sz);
	memcpy(&(this->rotation[0]),&(tmpRotation[0]),sz);
	memcpy(&(this->rotation[1]),&(tmpRotation[1]),sz);
	memcpy(&(this->rotation[2]),&(tmpRotation[2]),sz);
} 
/*****************************************************************************/
void Frame::mul(const Frame &f) 
{
	int i,j;
    double tmpRotation[3][3];
    double tmpTranslation[3];

    for (i=0;i<3;i++) {
			for (j=0;j<3;j++) {
				tmpRotation[i][j]=this->rotation[i][0]*f.rotation[0][j] +
					                this->rotation[i][1]*f.rotation[1][j] +
						              this->rotation[i][2]*f.rotation[2][j];
			}
    }
    for(i=0; i<3; i++) 
		tmpTranslation[i] = this->rotation[i][0]*f.translation[0] +
	                      this->rotation[i][1]*f.translation[1] +
	                      this->rotation[i][2]*f.translation[2] +
	                      this->translation[i];

	int sz = 3*sizeof(double);
	memcpy(this->translation,tmpTranslation,sz);
	memcpy(&(this->rotation[0]),&(tmpRotation[0]),sz);
	memcpy(&(this->rotation[1]),&(tmpRotation[1]),sz);
	memcpy(&(this->rotation[2]),&(tmpRotation[2]),sz);
}
/*****************************************************************************/
void Frame::invert() 
{
	double tmpTranslation[3];
	double tmpRotation[3][3];

	for(int i=0; i<3; i++)
	    for(int j=0; j<3; j++)
		tmpRotation[i][j] = this->rotation[j][i];

	tmpTranslation[0] = -(tmpRotation[0][0] * this->translation[0] + 
	                      tmpRotation[0][1] * this->translation[1] + 
	                      tmpRotation[0][2] * this->translation[2]);
	tmpTranslation[1] = -(tmpRotation[1][0] * this->translation[0] + 
	                      tmpRotation[1][1] * this->translation[1] + 
	                      tmpRotation[1][2] * this->translation[2]);
	tmpTranslation[2] = -(tmpRotation[2][0] * this->translation[0] + 
	                      tmpRotation[2][1] * this->translation[1] + 
	                      tmpRotation[2][2] * this->translation[2]);

	int sz = 3*sizeof(double);
	memcpy(this->translation,tmpTranslation,sz);
	memcpy(&(this->rotation[0]),&(tmpRotation[0]),sz);
	memcpy(&(this->rotation[1]),&(tmpRotation[1]),sz);
	memcpy(&(this->rotation[2]),&(tmpRotation[2]),sz);
}
/*****************************************************************************/
void Frame::invert(const Frame &f) 
{
	for(int i=0; i<3; i++)
	  for(int j=0; j<3; j++)
			this->rotation[i][j] = f.rotation[j][i];
	this->translation[0] = -(f.rotation[0][0] * f.translation[0] + 
			                     f.rotation[1][0] * f.translation[1] + 
			                     f.rotation[2][0] * f.translation[2]);
	this->translation[1] = -(f.rotation[0][1] * f.translation[0] + 
			                     f.rotation[1][1] * f.translation[1] + 
			                     f.rotation[2][1] * f.translation[2]);
	this->translation[2] = -(f.rotation[0][2] * f.translation[0] + 
			                     f.rotation[1][2] * f.translation[1] + 
			                     f.rotation[2][2] * f.translation[2]);
}
/*****************************************************************************/
void Frame::lerp(const Frame &f, double t, Frame &result) 
{
	double thisQuat[4], thisTrans[3], otherQuat[4], otherTrans[3];
	double lerpQuat[4], lerpTrans[3];

	getTranslation(thisTrans);
	getRotationQuaternion(thisQuat);
	f.getTranslation(otherTrans);
	f.getRotationQuaternion(otherQuat);

	lerpQuat[0] = (1-t) * thisQuat[0] + t * otherQuat[0];
	lerpQuat[1] = (1-t) * thisQuat[1] + t * otherQuat[1];
	lerpQuat[2] = (1-t) * thisQuat[2] + t * otherQuat[2];
	lerpQuat[3] = (1-t) * thisQuat[3] + t * otherQuat[3];

	lerpTrans[0] = (1-t) * thisTrans[0] + t * otherTrans[0];
	lerpTrans[1] = (1-t) * thisTrans[1] + t * otherTrans[1];
	lerpTrans[2] = (1-t) * thisTrans[2] + t * otherTrans[2];

	result.setRotationQuaternion(lerpQuat,true);
	result.setTranslation(lerpTrans);
}
/*****************************************************************************/
void Frame::lerp(const Frame &f, std::vector<double> &t, std::vector<Frame> &result)
{
	double thisQuat[4], thisTrans[3], otherQuat[4], otherTrans[3];
	double lerpQuat[4], lerpTrans[3];
	Frame res;

	result.clear();
	getTranslation(thisTrans);
	getRotationQuaternion(thisQuat);
	f.getTranslation(otherTrans);
	f.getRotationQuaternion(otherQuat);

	int sz = static_cast<int>(t.size());
	for(int i=0; i<sz; i++) {
		lerpQuat[0] = (1-t[i]) * thisQuat[0] + t[i] * otherQuat[0];
		lerpQuat[1] = (1-t[i]) * thisQuat[1] + t[i] * otherQuat[1];
		lerpQuat[2] = (1-t[i]) * thisQuat[2] + t[i] * otherQuat[2];
		lerpQuat[3] = (1-t[i]) * thisQuat[3] + t[i] * otherQuat[3];
		
		lerpTrans[0] = (1-t[i]) * thisTrans[0] + t[i] * otherTrans[0];
		lerpTrans[1] = (1-t[i]) * thisTrans[1] + t[i] * otherTrans[1];
		lerpTrans[2] = (1-t[i]) * thisTrans[2] + t[i] * otherTrans[2];
		
		res.setRotationQuaternion(lerpQuat,true);
		res.setTranslation(lerpTrans);
		result.push_back(res);
	}
}
/*****************************************************************************/
void Frame::slerp(const Frame &f, double t, Frame &result) 
{
	double thisQuat[4], thisTrans[3], otherQuat[4], otherTrans[3];
	double slerpQuat[4], slerpTrans[3];
	double dotProd, theta,thisScaleFactor, otherScaleFactor;

	getTranslation(thisTrans);
	getRotationQuaternion(thisQuat);
	f.getTranslation(otherTrans);
	f.getRotationQuaternion(otherQuat);
	
	dotProd = thisQuat[0]*otherQuat[0] + 
			  thisQuat[1]*otherQuat[1] +
			  thisQuat[2]*otherQuat[2] +
			  thisQuat[3]*otherQuat[3];
	theta = acos(dotProd);


	thisScaleFactor = sin((1-t)*theta) / sin(theta);
	otherScaleFactor = sin(t*theta) / sin(theta);

  slerpQuat[0] = thisScaleFactor * thisQuat[0] + otherScaleFactor * otherQuat[0];
  slerpQuat[1] = thisScaleFactor * thisQuat[1] + otherScaleFactor * otherQuat[1];
  slerpQuat[2] = thisScaleFactor * thisQuat[2] + otherScaleFactor * otherQuat[2];
	slerpQuat[3] = thisScaleFactor * thisQuat[3] + otherScaleFactor * otherQuat[3];

	slerpTrans[0] = (1-t) * thisTrans[0] + t * otherTrans[0];
	slerpTrans[1] = (1-t) * thisTrans[1] + t * otherTrans[1];
	slerpTrans[2] = (1-t) * thisTrans[2] + t * otherTrans[2];

	result.setRotationQuaternion(slerpQuat,false);
	result.setTranslation(slerpTrans);
}
/*****************************************************************************/
void Frame::slerp(const Frame &f, std::vector<double> &t, std::vector<Frame> &result)
{
	double thisQuat[4], thisTrans[3], otherQuat[4], otherTrans[3];
	double slerpQuat[4], slerpTrans[3];
	Frame res;
	double dotProd, theta,thisScaleFactor, otherScaleFactor, scaleFactor;

	result.clear();
	getTranslation(thisTrans);
	getRotationQuaternion(thisQuat);
	f.getTranslation(otherTrans);
	f.getRotationQuaternion(otherQuat);

	dotProd = thisQuat[0]*otherQuat[0] + 
			      thisQuat[1]*otherQuat[1] +
			      thisQuat[2]*otherQuat[2] +
			      thisQuat[3]*otherQuat[3];
	theta = acos(dotProd);

	scaleFactor = 1.0/sin(theta);

	int sz = static_cast<int>(t.size());
	for(int i=0; i<sz; i++) {
		thisScaleFactor = sin((1-t[i])*theta) * scaleFactor;
		otherScaleFactor = sin(t[i]*theta) * scaleFactor;

	  slerpQuat[0] = thisScaleFactor * thisQuat[0] + otherScaleFactor * otherQuat[0];
		slerpQuat[1] = thisScaleFactor * thisQuat[1] + otherScaleFactor * otherQuat[1];
		slerpQuat[2] = thisScaleFactor * thisQuat[2] + otherScaleFactor * otherQuat[2];
		slerpQuat[3] = thisScaleFactor * thisQuat[3] + otherScaleFactor * otherQuat[3];
		
		slerpTrans[0] = (1-t[i]) * thisTrans[0] + t[i] * otherTrans[0];
		slerpTrans[1] = (1-t[i]) * thisTrans[1] + t[i] * otherTrans[1];
		slerpTrans[2] = (1-t[i]) * thisTrans[2] + t[i] * otherTrans[2];
		
		res.setRotationQuaternion(slerpQuat,false);
		res.setTranslation(slerpTrans);
		result.push_back(res);
	}
}
/*****************************************************************************/
void Frame::setIdentity() 
{
	int sz = 3*sizeof(double);
	memset(this->translation,'\0',sz);
	memset(&(this->rotation[0]),'\0',sz);
	memset(&(this->rotation[1]),'\0',sz);
	memset(&(this->rotation[2]),'\0',sz);
	this->rotation[0][0] = 1;
	this->rotation[1][1] = 1;
	this->rotation[2][2] = 1;
}
/*****************************************************************************/
void Frame::setTranslation(double x, double y, double z) 
{
	this->translation[0] = x;
	this->translation[1] = y;
	this->translation[2] = z;
}
/*****************************************************************************/
void Frame::setTranslation(double translation[3]) 
{
	memcpy(this->translation,translation,3*sizeof(double));
}
/*****************************************************************************/
#ifdef USING_VNL
void Frame::setTranslation(vnl_vector<double> &translation) 
{
	this->translation[0] = translation[0];
	this->translation[1] = translation[1];
	this->translation[2] = translation[2];	
}
#endif //USING_VNL
/*****************************************************************************/
void Frame::setRotationAngles(double ax, double ay, double az) 
{
	double cx, cy, cz, sx, sy, sz;

	cx = cos(ax);
	cy = cos(ay);
	cz = cos(az);
	sx = sin(ax);
	sy = sin(ay);
	sz = sin(az);

	this->rotation[0][0] = cz*cy; 
	this->rotation[0][1] = cz*sy*sx - sz*cx;  
	this->rotation[0][2] = cz*sy*cx+sz*sx;     

	this->rotation[1][0] = sz*cy; 
	this->rotation[1][1] = sz*sy*sx + cz*cx ; 
	this->rotation[1][2] = sz*sy*cx - cz*sx;   

	this->rotation[2][0] = -sy;   
	this->rotation[2][1] = cy*sx;             
	this->rotation[2][2] = cy*cx;
}
/*****************************************************************************/
void Frame::setRotationAngles(double angles[3]) 
{
	double cx, cy, cz, sx, sy, sz;

	cx = cos(angles[0]);
	cy = cos(angles[1]);
	cz = cos(angles[2]);
	sx = sin(angles[0]);
	sy = sin(angles[1]);
	sz = sin(angles[2]);

	this->rotation[0][0] = cz*cy; 
	this->rotation[0][1] = cz*sy*sx - sz*cx;  
	this->rotation[0][2] = cz*sy*cx+sz*sx;     

	this->rotation[1][0] = sz*cy; 
	this->rotation[1][1] = sz*sy*sx + cz*cx ; 
	this->rotation[1][2] = sz*sy*cx - cz*sx;   

	this->rotation[2][0] = -sy;   
	this->rotation[2][1] = cy*sx;             
	this->rotation[2][2] = cy*cx;
}
/*****************************************************************************/
void Frame::setRotationMatrix(double rotationMatrix[3][3])
{
	int sz = 3*sizeof(double);
	memcpy(&(this->rotation[0]),&(rotationMatrix[0]),sz);
	memcpy(&(this->rotation[1]),&(rotationMatrix[1]),sz);
	memcpy(&(this->rotation[2]),&(rotationMatrix[2]),sz);
}
/*****************************************************************************/
#ifdef USING_VNL
void Frame::setRotationMatrix(vnl_matrix<double> &rotationMatrix)
{
	this->rotation[0][0] = rotationMatrix[0][0];   this->rotation[0][1] = rotationMatrix[0][1];   this->rotation[0][2] = rotationMatrix[0][2];
	this->rotation[1][0] = rotationMatrix[1][0];   this->rotation[1][1] = rotationMatrix[1][1];   this->rotation[1][2] = rotationMatrix[1][2];
	this->rotation[2][0] = rotationMatrix[2][0];   this->rotation[2][1] = rotationMatrix[2][1];   this->rotation[2][2] = rotationMatrix[2][2];	
}
#endif //USING_VNL
/*****************************************************************************/
void Frame::setRotationMatrix(double m00, double m01, double m02, 
		                          double m10, double m11, double m12, 
					                    double m20, double m21, double m22)
{
	this->rotation[0][0] = m00;   this->rotation[0][1] = m01;   this->rotation[0][2] = m02;
	this->rotation[1][0] = m10;   this->rotation[1][1] = m11;   this->rotation[1][2] = m12;
	this->rotation[2][0] = m20;   this->rotation[2][1] = m21;   this->rotation[2][2] = m22;
}
/*****************************************************************************/
void Frame::setRotationAxisAngle(double axis[3], double angle)
{
	double c,s,v;

	c =  cos(angle);
	s =  sin(angle);
	v = 1-c;

	this->rotation[0][0] = axis[0]*axis[0]*v + c; 
	this->rotation[0][1] = axis[0]*axis[1]*v - axis[2]*s;  
	this->rotation[0][2] = axis[0]*axis[2]*v + axis[1]*s;     

	this->rotation[1][0] = axis[0]*axis[1]*v + axis[2]*s; 
	this->rotation[1][1] = axis[1]*axis[1]*v + c; 
	this->rotation[1][2] = axis[1]*axis[2]*v - axis[0]*s;   
	
	this->rotation[2][0] = axis[0]*axis[2]*v - axis[1]*s;   
	this->rotation[2][1] = axis[1]*axis[2]*v + axis[0]*s;             
	this->rotation[2][2] = axis[2]*axis[2]*v + c;
}
/*****************************************************************************/
void Frame::setRotationAxisAngle(double axisAngle[3]) 
{
	double c,s,v,angle;
	double axis[3];

	angle = sqrt(axisAngle[0]*axisAngle[0] +
         		 axisAngle[1]*axisAngle[1] +
		         axisAngle[2]*axisAngle[2]);
	axis[0] = axisAngle[0]/angle;
	axis[1] = axisAngle[1]/angle;
	axis[2] = axisAngle[2]/angle;

	c =  cos(angle);
	s =  sin(angle);
	v = 1-c;

	this->rotation[0][0] = axis[0]*axis[0]*v + c; 
	this->rotation[0][1] = axis[0]*axis[1]*v - axis[2]*s;  
	this->rotation[0][2] = axis[0]*axis[2]*v + axis[1]*s;     

	this->rotation[1][0] = axis[0]*axis[1]*v + axis[2]*s; 
	this->rotation[1][1] = axis[1]*axis[1]*v + c; 
	this->rotation[1][2] = axis[1]*axis[2]*v - axis[0]*s;   
	
	this->rotation[2][0] = axis[0]*axis[2]*v - axis[1]*s;   
	this->rotation[2][1] = axis[1]*axis[2]*v + axis[0]*s;             
	this->rotation[2][2] = axis[2]*axis[2]*v + c;
}
/*****************************************************************************/
void Frame::setRotationQuaternion(double s, double qx, double qy, double qz, bool normalizeQuaternion) 
{
	if(normalizeQuaternion) {
		double norm = sqrt(s*s + qx*qx + qy*qy + qz*qz);
		s/=norm;
		qx/=norm;
		qy/=norm;
		qz/=norm;
	}

	this->rotation[0][0] = 1-2*(qy*qy+qz*qz);   
	this->rotation[0][1] = 2*(qx*qy-s*qz);    
	this->rotation[0][2] = 2*(qx*qz+s*qy);    
	
	this->rotation[1][0] = 2*(qx*qy+s*qz);     
	this->rotation[1][1] = 1-2*(qx*qx+qz*qz);  
	this->rotation[1][2] = 2*(qy*qz-s*qx);    
	
	this->rotation[2][0] = 2*(qx*qz-s*qy);     
	this->rotation[2][1] = 2*(qy*qz+s*qx);    
	this->rotation[2][2] = 1-2*(qx*qx+qy*qy);  
}
/*****************************************************************************/
void Frame::setRotationQuaternion(double quaternion[4], bool normalizeQuaternion) 
{
	double s,qx,qy,qz;
	s = quaternion[0];
	qx = quaternion[1];
	qy = quaternion[2];
	qz = quaternion[3];

	if(normalizeQuaternion) {
		double norm = sqrt(s*s + qx*qx + qy*qy + qz*qz);
		s/=norm;
		qx/=norm;
		qy/=norm;
		qz/=norm;
	}

	this->rotation[0][0] = 1-2*(qy*qy+qz*qz);   
	this->rotation[0][1] = 2*(qx*qy-s*qz);    
	this->rotation[0][2] = 2*(qx*qz+s*qy);    
	
	this->rotation[1][0] = 2*(qx*qy+s*qz);     
	this->rotation[1][1] = 1-2*(qx*qx+qz*qz);  
	this->rotation[1][2] = 2*(qy*qz-s*qx);    
	
	this->rotation[2][0] = 2*(qx*qz-s*qy);     
	this->rotation[2][1] = 2*(qy*qz+s*qx);    
	this->rotation[2][2] = 1-2*(qx*qx+qy*qy);  
}
/*****************************************************************************/
void Frame::set(const Frame &f) 
{
	int sz = 3*sizeof(double);
	memcpy(this->translation,f.translation,sz);
	memcpy(&(this->rotation[0]),&(f.rotation[0]),sz);
	memcpy(&(this->rotation[1]),&(f.rotation[1]),sz);
	memcpy(&(this->rotation[2]),&(f.rotation[2]),sz);
}
/*****************************************************************************/
void Frame::getTranslation(double translation[3]) const
{
	memcpy(translation,this->translation,3*sizeof(double));
}
/*****************************************************************************/
void Frame::getTranslation(double &x, double &y, double &z) const
{
	x = this->translation[0];
	y = this->translation[1];
	z = this->translation[2];
}
/*****************************************************************************/
#ifdef USING_VNL
void Frame::getTranslation(vnl_vector<double> &translation) const
{
	translation[0] = this->translation[0];
	translation[1] = this->translation[1];
	translation[2] = this->translation[2];
}
#endif //USING_VNL
/*****************************************************************************/
void Frame::getRotationAngles(double angles[6], bool &isGimbalLock) const
{
  double cy;

	angles[1] = atan2(-this->rotation[2][0],
			              sqrt(this->rotation[0][0]*this->rotation[0][0] +
				                 this->rotation[1][0]*this->rotation[1][0]));
  angles[4] = atan2(-this->rotation[2][0],
			              -sqrt(this->rotation[0][0]*this->rotation[0][0] +
				                 this->rotation[1][0]*this->rotation[1][0]));

	if(fabs(angles[1] - Frame::HALF_PI) > Frame::SMALL_ANGLE &&
	   fabs(angles[1] + Frame::HALF_PI) > Frame::SMALL_ANGLE) {
		isGimbalLock = false;
    cy = cos(angles[1]);
	  angles[0] = atan2(this->rotation[2][1]/cy, 
					            this->rotation[2][2]/cy);
	  angles[2] = atan2(this->rotation[1][0]/cy, 
					            this->rotation[0][0]/cy);
	  cy = cos(angles[4]);
	  angles[3] = atan2(this->rotation[2][1]/cy, 
					            this->rotation[2][2]/cy);
	  angles[5] = atan2(this->rotation[1][0]/cy, 
					            this->rotation[0][0]/cy);
	}
          	//thetaY is approximatly PI or -PI
           //set thetaZ to zero and compute thetaX
	else {
		isGimbalLock = true;
	    angles[2] = angles[5] = 0;
	    angles[0] = angles[3] = atan2(this->rotation[0][1],this->rotation[1][1]);
	}

}
/*****************************************************************************/
void Frame::getRotationAxisAngle(double axisAngle[4]) const 
{
	double cTheta = (this->rotation[0][0] + this->rotation[1][1] + this->rotation[2][2] - 1) / 2.0;
	double sTheta = sqrt(((this->rotation[1][0] - this->rotation[0][1])*(this->rotation[1][0] - this->rotation[0][1]) +
		                   (this->rotation[0][2] - this->rotation[2][0])*(this->rotation[0][2] - this->rotation[2][0]) +
		                   (this->rotation[2][1] - this->rotation[1][2])*(this->rotation[2][1] - this->rotation[1][2]))/4.0);

	axisAngle[0] = atan2(sTheta,cTheta);
	          //not a rotation angle near zero or pi
	if((axisAngle[0] > Frame::SMALL_ANGLE) &&
	   (axisAngle[0] < (Frame::PI - Frame::SMALL_ANGLE))) {
		double scale = (1.0/(2*sin(axisAngle[0])));
	    axisAngle[1] = scale*(this->rotation[2][1] - this->rotation[1][2]);
	    axisAngle[2] = scale*(this->rotation[0][2] - this->rotation[2][0]);
	    axisAngle[3] = scale*(this->rotation[1][0] - this->rotation[0][1]);	
	}   //stablize the computation of the axis when sin(axisAngle[0]) is close to zero 
	else {  
		int i,j,k;
		double w;
		       //find maximal entry on diagonal
		i=0;
		if(this->rotation[1][1] > this->rotation[i][i])
			i=1;
		if(this->rotation[2][2] > this->rotation[i][i])
			i=2;
		j = (i+1)%3;
		k = (j+1)%3;

		w = 1.0/(2*(1 - cTheta));
		       //compute the vector part
		axisAngle[i+1] = sqrt((this->rotation[i][i] - this->rotation[j][j] - this->rotation[k][k] + 1)*w);
		axisAngle[j+1] = (this->rotation[i][j] + this->rotation[j][i])*(w/axisAngle[i+1]);
		axisAngle[k+1] = (this->rotation[i][k] + this->rotation[k][i])*(w/axisAngle[i+1]);
	}
}
/*****************************************************************************/
void Frame::getRotationQuaternion(double quaternion[4]) const
{
	double startSingularRange = Frame::HALF_PI - Frame::SMALL_ANGLE;
	double endSingularRange = Frame::HALF_PI + Frame::SMALL_ANGLE;

	quaternion[0] = (0.5*sqrt(this->rotation[0][0] +
					                  this->rotation[1][1] +
					                  this->rotation[2][2] +
					                  1));
	          //when s is near zero then halfTheta is near 90 degrees
	double halfTheta = acos(quaternion[0]);
	          //check that we are not in the singular zone
	if(!(halfTheta > startSingularRange   && halfTheta < endSingularRange)) {
	    double denom = 4*quaternion[0];
	    quaternion[1] = (this->rotation[2][1] - this->rotation[1][2])/denom;
	    quaternion[2] = (this->rotation[0][2] - this->rotation[2][0])/denom;
	    quaternion[3] = (this->rotation[1][0] - this->rotation[0][1])/denom;
	}
	else {  //stablize the computation of the vector part when halfTheta is near 90 degrees
		int i,j,k;
		double w;
		       //find maximal entry on diagonal
		i=0;
		if(this->rotation[1][1] > this->rotation[i][i])
			i=1;
		if(this->rotation[2][2] > this->rotation[i][i])
			i=2;
		j = (i+1)%3;
		k = (j+1)%3;

		       //compute the vector part
		w = sqrt(this->rotation[i][i] - this->rotation[j][j] - this->rotation[k][k] + 1);
		quaternion[i+1] = w/2.0;
		quaternion[j+1] = (this->rotation[i][j] + this->rotation[j][i])/(2*w);
		quaternion[k+1] = (this->rotation[i][k] + this->rotation[k][i])/(2*w);
	}
}
/*****************************************************************************/
void Frame::getRotationMatrix(double rotationMatrix[3][3]) const
{
	int sz = 3*sizeof(double);
	memcpy(&(rotationMatrix[0]),&(this->rotation[0]),sz);
	memcpy(&(rotationMatrix[1]),&(this->rotation[1]),sz);
	memcpy(&(rotationMatrix[2]),&(this->rotation[2]),sz);
}
/*****************************************************************************/
void Frame::getRotationMatrix(double &m00, double &m01, double &m02, 
		                          double &m10, double &m11, double &m12, 
					                    double &m20, double &m21, double &m22) const
{
	m00 = this->rotation[0][0];   m01 = this->rotation[0][1];   m02 = this->rotation[0][2];
	m10 = this->rotation[1][0];   m11 = this->rotation[1][1];   m12 = this->rotation[1][2];
	m20 = this->rotation[2][0];   m21 = this->rotation[2][1];   m22 = this->rotation[2][2];
}
/*****************************************************************************/
#ifdef USING_VNL
void Frame::getRotationMatrix(vnl_matrix<double> &rotationMatrix)
{
	rotationMatrix[0][0] = this->rotation[0][0];   rotationMatrix[0][1] = this->rotation[0][1];   rotationMatrix[0][2] = this->rotation[0][2];
	rotationMatrix[1][0] = this->rotation[1][0];   rotationMatrix[1][1] = this->rotation[1][1];   rotationMatrix[1][2] = this->rotation[1][2];
	rotationMatrix[2][0] = this->rotation[2][0];   rotationMatrix[2][1] = this->rotation[2][1];   rotationMatrix[2][2] = this->rotation[2][2];	
}
#endif //USING_VNL
/*****************************************************************************/
bool Frame::angleAndTranslationDiff(const Frame &f, double &dx, double &dy, double &dz, 
									                  double &dax, double &day, double &daz) const
{
	double thisRotationAngles[6], thisTranslation[3], otherRotationAngles[6], otherTranslation[3];
	bool isGimbalLock;

	getTranslation(thisTranslation);
	getRotationAngles(thisRotationAngles, isGimbalLock);
	if(isGimbalLock)
		return false;

	f.getTranslation(otherTranslation);
	f.getRotationAngles(otherRotationAngles, isGimbalLock);
	if(isGimbalLock)
		return false;

	dx = fabs(thisTranslation[0] - otherTranslation[0]);
	dy = fabs(thisTranslation[1] - otherTranslation[1]);
	dz = fabs(thisTranslation[2] - otherTranslation[2]);
	dax = fabs(thisRotationAngles[0] - otherRotationAngles[0]);
	day = fabs(thisRotationAngles[1] - otherRotationAngles[1]);
	daz = fabs(thisRotationAngles[2] - otherRotationAngles[2]);
	return true;
}
/*****************************************************************************/
void Frame::angleAndTranslationDiff(const Frame &f, double &dx, double &dy, double &dz, 
							                      double &angle) const
{
	double axisAndAngle[4];
	Frame res;
	Frame fInverted;
	fInverted.invert(f);

	res.mul(fInverted,*this);
	
	res.getRotationAxisAngle(axisAndAngle);
	angle = fabs(axisAndAngle[0]);
	
	dx = fabs(res.translation[0]);
	dy = fabs(res.translation[1]);
	dz = fabs(res.translation[2]);	
}

} //namespace lsqrRecipes
