#ifndef _POINT_H_
#define _POINT_H_

#include "copyright.h"

#include <vnl/vnl_vector_ref.h> 
#include "Vector.h"

/**
 * Template defining nD points. Valid types for point entries are the numeric types that 
 * can be reasonably cast to double (int, float, etc.).
 *
 * @author: Ziv Yaniv (zivy@isis.georgetown.edu)
 */

namespace lsqrRecipes {

template <class T, unsigned int n>
class Point {
public:
              //enable compile time access to the point dimension, ugly but necessary 
	enum {dimension=n};

	/**
	 * Ouput the point data to the given output stream as the following string:
	 *[x1,x2,...xn]
	 */
	friend std::ostream &operator<<(std::ostream& output, const Point &p) {
		output<<"[ "<<p.data[0];
		if(n==1) 
			output<<" ]"; 
		else {
			for(unsigned int i=1; i<n; i++) 
				output<<", "<<p.data[i]; 
			output<<" ]";
		}
		return output;
	}

  /**
   * Assignment operator.
   */
  Point<T,n> &operator =(const Point<T,n> &other) {memcpy(data,other.data,n*sizeof(T)); return *this;}

	/**
	 * Access to the coordinates of this point. No bounds checking is performed.
   * Valid indexes are in [0,n-1].
	 */
	double & operator[](int index) {return this->data[index];}
  const double & operator[](int index) const {return this->data[index];}

  /**
	 * Default constructor, sets the point to zero.
	 */
	Point(){memset(this->data,0,n*sizeof(T));}

  /**
	 * Construct a point from the given pointer to an array.
	 * @param fillData An array that has at least 'n' elements.
	 */
	Point(T *fillData) {memcpy(this->data,fillData,n*sizeof(T));}

  /**
	 * Copy constructor.
	 * @param other The point we copy.
	 */
	Point(const Point<T,n> &other) {memcpy(this->data,other.data,n*sizeof(T));}
			
  /**
   * Return the distance between this and the given point.
   */
  double distanceSquared(const Point<T,n> &other) {
    return (this->getVnlVector() - other.getVnlVector()).squared_magnitude();
  }

  /**
	 * Fill the point from the given pointer to an array.
	 * @param fillData An array that has at least 'n' elements.
	 */
	void set(T *fillData) {memcpy(this->data,fillData,n*sizeof(T));}

	/**
	 * Return the dimensionality of this point.
	 */
	unsigned int size() {return n;}

  /**
   * Add a vector to the current point.
   * @param v The point is translated using this vector.
   */
  Point operator+(const Vector<T,n> &v) {
    Point result(*this);
    for(unsigned int i=0; i<n; i++)
      result[i]+=v[i];
    return result;
  }

  /**
   * Subtract a vector from the current point.
   * @param v The point is translated using this vector.
   */
  Point operator-(const Vector<T,n> &v) {
    Point result(*this);
    for(unsigned int i=0; i<n; i++)
      result[i]-=v[i];
    return result;
  }

  /**
   * Subtract points to get vector.
   * @param p The point to subtract.
   */
  Vector<T,n> operator-(const Point &p) {
    Vector<T,n> result(this->data);
    for(unsigned int i=0; i<n; i++)
      result[i]-=p[i];
    return result;
  }

private:

	/**
	 * Return this point as a vnl_vector_ref object.
	 */
	vnl_vector_ref<T> getVnlVector() const {return vnl_vector_ref<T>(n, const_cast<double *>(this->data));}

  T data[n];
};

} //namespace lsqrRecipes

#endif //_POINT_H_
