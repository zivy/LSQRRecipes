#ifndef _FRAME_H_
#define _FRAME_H_

#include "copyright.h"

#include <vector>
#include "Point3D.h"
#include "Vector3D.h"

#define USING_VNL

#ifdef USING_VNL
	#include <vnl/vnl_matrix.h>
	#include <vnl/vnl_vector.h>
  #include <vnl/vnl_vector_fixed.h> 
#endif

/**
 * Implementation of a Frame/Rigid/Euclidean transformation.
 * All angular input and ouput is done in radians.
 *
 * @author: Ziv Yaniv (zivy@isis.georgetown.edu)
 */

namespace lsqrRecipes {

class Frame {
	friend std::ostream &operator<<(std::ostream& output, const Frame &f);
private:
	double rotation[3][3];
	double translation[3]; 
public:                                    //0.5 degrees
  static const double SMALL_ANGLE;
  static const double PI;
  static const double HALF_PI;
	static const double DEGREES_PER_RADIAN;
	static const double RADIANS_PER_DEGREE;
	

	enum {MATRIX=0, EULER_ANGLES, AXIS_ANGLE, QUATERNION};
  int outputFormat;

	/**
	 * Usefull method when setting frame data. 
	 * Translates degrees to radians.
	 */
	static inline double toRadians(double angleInDegrees) {
		return angleInDegrees * RADIANS_PER_DEGREE;
	}

	/**
	 * Usefull method when getting frame data. 
	 * Translates radians to degrees.
	 */
	static inline double toDegrees(double angleInRadians) {
		return angleInRadians * DEGREES_PER_RADIAN;
	}


    /**
     * Create a frame with the given rotation matrix and translation, if
		 * both parameters are NULL then create the identity transformation.
     * @param R Rotation matrix (3X3) in row major order (R[i][j], i'th row j'th column).
     * @param t Translation vector (3X1).
		 */
	Frame(double R[3][3] = NULL, double t[3] = NULL);


	  /**
     * Copy Constructor, create a frame which is a copy of the given frame.
     * @param f This frame is copied.
     */
	Frame(const Frame &f);


    /**
     * Create a frame with the given translations and rotation angles.
     * @param x Translation in the x direction.
     * @param y Translation in the y direction.
     * @param z Translation in the z direction.
     * @param az Angle of rotation around the z axis.
     * @param ay Angle of rotation around the y axis.
     * @param ax Angle of rotation around the x axis.
     */
	Frame(double x, double y, double z, double ax, double ay, double az);


    /**
     * Create a frame with the given translations and rotation around the
     * given axis with the given angle. Axis of rotation is assumed to
     * be normalized.
     * @param x Translation in the x direction.
     * @param y Translation in the y direction.
     * @param z Translation in the z direction.
     * @param axis Axis of rotation.
     * @param angle Angle of rotation around the given axis.
     */
	Frame(double x, double y, double z, double axis[3], double angle);


    /**
     * Create a frame with the given translations and rotation around the
     * given axis and angle. The norm of the axis is the size of the
     * rotation angle.
     * @param x Translation in the x direction.
     * @param y Translation in the y direction.
     * @param z Translation in the z direction.
     * @param axisAngle The direction of this vector is the axis of rotation
     *                  and its norm is the size of the rotation angle.
     */
	Frame(double x, double y, double z, double axisAngle[3]);


    /**
     * Create a frame with the given translations and rotation according
     * to the given quaternion.
     * @param x Translation in the x direction.
     * @param y Translation in the y direction.
     * @param z Translation in the z direction.
     * @param s Scalar part of the quaternion.
     * @param qx x coordinate of the quaternion.
     * @param qy y coordinate of the quaternion.
     * @param qz z coordinate of the quaternion.
	   * @param normalizeQuaternion true if we need to normalize the quaternion.
	   *                            Rotations are unit quaternions (s^2 + qx^2 + qy^2 + qz^2 = 1).
     */
	Frame(double x, double y, double z, double s, double qx, double qy, double qz, bool normalizeQuaternion = false);

	/**
	 * Set the output format to the given parameter and return the previous output format.
	 * @param newFormat The format in which the frame will write itself to output streams.
	 * @return Returns the previous format.
	 */
	int setOutputFormat(int newFormat);

    /**
     * Apply the transformation on the given point and update the
     * point coordinates accordingly.
     * @param p The point on which the transformation is applied.
     * @see Point3D
     */
	void apply(Point3D &p) const;


    /**
     * Apply the transformation on the given point and set the results in
     * the second point.
     * @param p The point on which the transformation is applied.
     * @param pTransformed The result of applying the transformation on point
     *                     p.
     * @see Point3D
     */
	void apply(const Point3D &p, Point3D &pTransformed) const;


    /**
     * Apply the inverse transformation on the given point and update the
     * point coordinates accordingly.
     * If the inverse transformation is applied many times 
     * it is computationally more efficient to invert the frame and 
     * then use the inverted frame's apply() method.
     * @param p The point on which the transformation is applied.
     * @see Point3D
     */
	void applyInverse(Point3D &p);


    /**
     * Apply the inverse transformation on the given point and update the
     * point coordinates accordingly.
     * If the inverse transformation is applied many times 
     * it is computationally more efficient to invert the frame and 
     * then use the inverted frame's apply() method.
     * @param p The point on which the transformation is applied.
     * @param pTransformed The result of applying the transformation on point
     *                     p.
     * @see Point3D
     */
	void applyInverse(Point3D &p, Point3D &pTransformed);


    /**
     * Apply the transformation on the given vector and update the
     * vector coordinates accordingly.
		 * The resulting vector is not affected by the translational part of the
		 * transformation, only by the rotational part.
     * @param v The vector on which the transformation is applied.
     * @see Vector3D
     */
	void apply(Vector3D &v) const;
#ifdef USING_VNL
  void apply(vnl_vector_fixed<double,3> &v) const;
#endif

    /**
     * Apply the transformation on the given vector and set the results in
     * the second vector.
		 * The resulting vector is not affected by the translational part of the
		 * transformation, only by the rotational part.
     * @param v The vector on which the transformation is applied.
     * @param vTransformed The result of applying the transformation on vector
     *                     v.
		 *
     * @see Vector3D
     */
	void apply(Vector3D &v, Vector3D &vTransformed) const;
#ifdef USING_VNL
  void apply(vnl_vector_fixed<double,3> &v, vnl_vector_fixed<double,3> &vTransformed) const;
#endif


    /**
     * Apply the inverse transformation on the given vector and update the
     * vector coordinates accordingly.
     * If the inverse transformation is applied many times 
     * it is computationally more efficient to invert the frame and 
     * then use the inverted frame's apply() method.
		 * The resulting vector is not affected by the translational part of the
		 * transformation, only by the rotational part.
		 * @param v The vector on which the transformation is applied.
     * @see Vector3D
     */
	void applyInverse(Vector3D &v);
#ifdef USING_VNL
  void applyInverse(vnl_vector_fixed<double,3> &v);
#endif


    /**
     * Apply the inverse transformation on the given vector and update the
     * vector coordinates accordingly.
     * If the inverse transformation is applied many times 
     * it is computationally more efficient to invert the frame and 
     * then use the inverted frame's apply() method.
		 * The resulting vector is not affected by the translational part of the
		 * transformation, only by the rotational part.
     * @param v The point on which the transformation is applied.
     * @param vTransformed The result of applying the transformation on point
     *                     v.
     * @see Point3D
     */
	void applyInverse(Vector3D &v, Vector3D &vTransformed);
#ifdef USING_VNL
	void applyInverse(vnl_vector_fixed<double,3> &v, vnl_vector_fixed<double,3> &vTransformed); 
#endif

  	/**
     * Compute the composition/multiplication (f1 o f2)
     * of two frames and store the result in this frame.
     * @param f1 Left operand of the composition.
     * @param f2 Right operand of the composition.
     */
  void mul(const Frame &f1, const Frame &f2);


    /**
     * Compute the composition/multiplication (this o f)
     * of two frames and store the result in this frame.
     * @param f Right operand of the composition.
     */
  void mul(const Frame &f);


    /**
     * Invert the frame.
     */
  void invert();


    /**
     * Invert the given frame and store the result in this frame.
     * @param f The frame whose inverse we set this frame to be.
     */
	void invert(const Frame &f);


	/**
	 * Compute a frame which is between this frame and the given frame
	 * at the given parameter value (t in [0,1]) using linear 
	 * interpolation (lerp).
	 * @param f The next frame.
	 * @param t The interpolation parameter. Assumes t is in [0,1].
	 * @param result The result of the interpolation [(1-t)*this + t*f]. 
	 */
	void lerp(const Frame &f, double t, Frame &result);

	
	/**
	 * Compute frames which are between this frame and the given frame
	 * at the given parameter values (t in [0,1]) using linear
	 * interpolation (lerp).
	 * @param f The next frame.
	 * @param t The vector of interpolation parameters. Assumes t is in [0,1].
   * @param result A vector containing the interpolated frames. 
	 *               The vector is cleared by the function and then
	 *               filled with the interpolating frames [(1-t)*this + t*f].
	 */
	void lerp(const Frame &f, std::vector<double> &t, std::vector<Frame> &result);

	
	/**
	 * Compute a frame which is between this frame and the given frame
	 * at the given parameter value (t in [0,1]) using spherical linear
	 * interpolation (slerp).
	 * @param f The next frame.
	 * @param t The interpolation parameter. Assumes t is in [0,1].
	 * @param result The result of the interpolation 
	 *        sin((1-t)*Theta)            sin(t*Theta)
	 *        ----------------- * this +  ------------- * f
	 *             sin(Theta)               sin(Theta)
	 */
	void slerp(const Frame &f, double t, Frame &result);

	
	/**
	 * Compute frames which are between this frame and the given frame
	 * at the given parameter values (t in [0,1]) using spherical linear
	 * interpolation (slerp).
	 * @param f The next frame.
	 * @param t The vector of interpolation parameters. Assumes t is in [0,1].
   * @param result A vector containing the interpolated frames. 
	 *               The vector is cleared by the function and then
	 *               filled with the interpolating frames.
	 *        sin((1-t)*Theta)            sin(t*Theta)
	 *        ----------------- * this +  ------------- * f
	 *             sin(Theta)               sin(Theta)
	 */
	void slerp(const Frame &f, std::vector<double> &t, std::vector<Frame> &result);

    
	  /**
     * Set this frame to the identity transformation, no rotation and no
     * translation.
     */
  void setIdentity();


    /**
     * Set this frame's translation to the given translation.
     * @param x Translation in the x direction.
     * @param y Translation in the y direction.
     * @param z Translation in the z direction.
     */
  void setTranslation(double x, double y, double z);


    /**
     * Set this frame's translation to the given translation.
     * Assumes that the given array has at least three entries.
     * @param translation Array containing the translation parameters
     *                    [x,y,z].
     */
  void setTranslation(double translation[3]);

#ifdef USING_VNL
	  /**
     * Set this frame's translation to the given translation.
     * Assumes that the given vnl_vector has at least three entries.
     * @param translation vnl_vector<double> containing the translation parameters
     *                    [x,y,z].
     */
  void setTranslation(vnl_vector<double> &translation);
#endif

    /**
     * Set this frame's rotation according to the given angles.
     * @param az Angle of rotation around the z axis.
     * @param ay Angle of rotation around the y axis.
     * @param ax Angle of rotation around the x axis.
     */
  void setRotationAngles(double ax, double ay, double az);


    /**
     * Set this frame's rotation according to the given angles.
     * @param angles Angles of rotation [ax,ay,az].
     */
  void setRotationAngles(double angles[3]);

    
	  /**
     * Set the rotation of this frame according to the given rotation matrix.
		 * Does not validate that the matrix is a rotation (i.e. RR^T = I and det(R) = 1)
     * @param rotationMatrix An array containing the rotation matrix parameters.
     *                           [m00 m01 m02].
	   *                           [m10 m11 m12]
	   *                           [m20 m21 m22]
     */
  void setRotationMatrix(double rotationMatrix[3][3]);

#ifdef USING_VNL
	  /**
     * Set the rotation of this frame according to the given rotation matrix.
		 * Does not validate that the matrix is a rotation (i.e. RR^T = I and det(R) = 1)
     * @param rotationMatrix A vnl_matrix<double> containing the rotation matrix parameters.
     *                           [m00 m01 m02].
	   *                           [m10 m11 m12]
	   *                           [m20 m21 m22]
     */
  void setRotationMatrix(vnl_matrix<double> &rotationMatrix);
#endif

    /**
     * Set the rotation of this frame according to the given rotation matrix.
		 * Does not validate that the matrix is a rotation (i.e. RR^T = I and det(R) = 1)
     * @param m(i,j) The entries of the matrix.
	   */
	void setRotationMatrix(double m00, double m01, double m02, 
		                     double m10, double m11, double m12, 
			     		           double m20, double m21, double m22);


    /**
     * Set this frame's rotation according to the given axis and angle.
     * Axis of rotation is assumed to be normalized.
     * @param axis Axis of rotation.
     * @param angle Angle of rotation around the given axis.
     */ 
  void setRotationAxisAngle(double axis[3], double angle);


    /**
     * Set this frame's rotation according to the given axis and angle.
     * The norm of the axis is the size of the rotation angle.
     * @param axisAngle The direction of this vector is the axis of rotation
     *                  and its norm is the size of the rotation angle.
     */
  void setRotationAxisAngle(double axisAngle[3]);


    /**
     * Set this frame's rotation according to the given quaternion.
     * @param s Scalar part of the quaternion.
     * @param qx x coordinate of the quaternion.
     * @param qy y coordinate of the quaternion.
     * @param qz z coordinate of the quaternion.
	   * @param normalizeQuaternion true if we need to normalize the quaternion.
	   *                            Rotations are unit quaternions (s^2 + qx^2 + qy^2 + qz^2 = 1).
     */
  void setRotationQuaternion(double s, double qx, double qy, double qz, bool normalizeQuaternion = false);


    /**
     * Set this frame's rotation according to the given quaternion.
     * @param quat Quaternion data in the format [s,qx,qy,qz].
	   * @param normalizeQuaternion true if we need to normalize the quaternion.
	   *                            Rotations are unit quaternions (s^2 + qx^2 + qy^2 + qz^2 = 1).
     */
  void setRotationQuaternion(double quaternion[4], bool normalizeQuaternion = false);


    /**
     * Set this frame to the given frame.
     * @param f Set the current frame to this frame.
     */
  void set(const Frame &f);


    /**
     * Get the translational part of this frame.
     * @param translation An array into which the translation is 
     *                    written [x,y,z].
     */
  void getTranslation(double translation[3]) const;


    /**
     * Get the translational part of this frame.
     * @param x Translation in the x direction.
     * @param y Translation in the y direction.
     * @param z Translation in the z direction.
     */
  void getTranslation(double &x, double &y, double &z) const;
    

#ifdef USING_VNL
    /**
     * Get the translational part of this frame.
     * @param translation A vnl_vector<double> into which the translation is 
     *                    written [x,y,z].
     */
	void getTranslation(vnl_vector<double> &translation) const;
#endif //USING_VNL

	  /**
     * Get the angles of rotation from this frame. Rotation angles are
     * for the ZYX Euler angle specification.
     * @param angles An array into which the angles are  
     *                    written [ax1,ay1,az1,ax2,ay2,az2].
     *               In general there are two sets of angles that yield the 
     *               same matrix. If we have gimbal lock then there are an infinite
     *               number of angles (choices for ax,az). We arbitrarily set
     *               az to zero and compute ax.
     *               
	   * @param isGimbalLock Set to true if 1 or 2 (thetaY is approximatly PI or -PI).
	   *                     1. fabs(ay-90)<SMALL_ANGLE 
	   *                     2. fabs(ay+90)<SMALL_ANGLE
	   *                     In this case we cannot seperate the rotations 
	   *                     into ax,az. We arbitrarily set az=0 and ax is set accordingly.
     */
  void getRotationAngles(double angles[6], bool &isGimbalLock) const;


    /**
     * Get the axis and angle of rotation from this frame.
     * @param axisAngle An array into which the angle and the axis are
     *                  written [a,nx,ny,nz].
     */
  void getRotationAxisAngle(double axisAngle[4]) const;


    /**
     * Get the unit quaternion describing the rotation from this frame.
     * @param quaternion An array into which the quaternion parameters are
     *                   written [s,qx,qy,qz].
     */
  void getRotationQuaternion(double quaternion[4]) const;

    /**
     * Get the matrix describing the rotation from this frame.
     * @param rotationMatrix An array into which the matrix parameters are
     *                   written [m00 m01 m02].
	   *                           [m10 m11 m12]
	   *                           [m20 m21 m22]
     */
  void getRotationMatrix(double rotationMatrix[3][3]) const;


    /**
     * Get the matrix describing the rotation from this frame.
     * @param m(i,j) The entries of the matrix.
	   */
	void getRotationMatrix(double &m00, double &m01, double &m02, 
		                     double &m10, double &m11, double &m12, 
					               double &m20, double &m21, double &m22) const;

#ifdef USING_VNL
	  /**
     * Get the matrix describing the rotation from this frame.
     * @param rotationMatrix A vnl_matrix<double> whose size is at least 3X3 
		 *                       into which the matrix parameters are written into the top left corner.
     *                           [m00 m01 m02].
	   *                           [m10 m11 m12]
	   *                           [m20 m21 m22]
     */
	void getRotationMatrix(vnl_matrix<double> &rotationMatrix);
#endif //USING_VNL

	 /**
	  * Get the absolute difference in translation and euler rotation angles between this frame and the
	  * given one. If there is gimbal-lock in one of the frames return false otherwise true. 
	  */
	bool angleAndTranslationDiff(const Frame &f, double &dx, double &dy, double &dz, double &dax, double &day, double &daz) const;

	 /**
	  * Get the difference in translation between this frame and the
	  * given one.
	  * Essentialy, this shows how different the frame T = (f^-1 * this) is from the identity frame.
		* The difference in translation is just the absolute value of the transaltions of T.
		* The difference in rotation is given as a single angular value, this is the angle
		* from the axis angle representation of the rotational part of T.
    */
	void angleAndTranslationDiff(const Frame &f, double &dx, double &dy, double &dz, double &angle) const;
};

} //namespace lsqrRecipes

#endif  //_FRAME_H_
