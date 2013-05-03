#ifndef _VECTOR_H_
#define _VECTOR_H_

#include<string.h>
#include <vnl/vnl_vector_ref.h> 
#include "copyright.h"


/**
 * Template defining nD vectors. Valid types for vector entries are the numeric 
 * types that can be reasonably cast to double (int, float, etc.).
 *
 * @author: Ziv Yaniv
 */

namespace lsqrRecipes {

template <class T, unsigned int n>
class Vector {
public:
              //enable compile time access to the point dimension, ugly but necessary 
	enum {dimension=n};

	/**
	 * Ouput the vector data to the given output stream as the following string:
	 *[x1,x2,...xn]
	 */
	friend std::ostream &operator<<(std::ostream& output, const Vector &v) {
		output<<"[ "<<v.data[0];
		if(n==1) 
			output<<" ]"; 
		else {
			for(unsigned int i=1; i<n; i++) 
				output<<", "<<v.data[i]; 
			output<<" ]";
		}
		return output;
	}

  /**
   * Assignment operator.
   */
  Vector<T,n> &operator =(const Vector<T,n> &other) {memcpy(data,other.data,n*sizeof(T)); return *this;}

	/**
	 * Access to the coordinates of this vector. No bounds checking is performed.
   * Valid indexes are in [0,n-1].
	 */
	double & operator[](int index) {return this->data[index];}
  const double & operator[](int index) const {return this->data[index];}

  /**
	 * Default constructor, sets the vector to zero.
	 */
	Vector(){memset(this->data,0,n*sizeof(T));}

  /**
	 * Construct a vector from the given pointer to an array.
	 * @param fillData An array that has at least 'n' elements.
	 */
	Vector(T *fillData) {memcpy(this->data,fillData,n*sizeof(T));}

  /**
	 * Copy constructor.
	 * @param other The vector we copy.
	 */
	Vector(const Vector<T,n> &other) {memcpy(this->data,other.data,n*sizeof(T));}
	
  /**
   * Multiply vector by scalar.
   * @param scalar The vector is scaled with this value.
   */
  Vector  operator*(const T & scalar) const 
  {
    Vector result(*this);
    for(unsigned int i=0; i<n; i++)
      result[i]*=scalar;
    return result;
  }


  /**
   * Vector dot product.
   * @param right The right side of the dot product.
   */
  T operator*(const Vector & right) const 
  {
    T result;
    result = 0;
    for(unsigned int i=0; i<n; i++)
      result+= data[i]*right[i];
    return result;
  }


  /**
   * Vector addition.
   * @param right The right operand of the plus sign.
   */
  Vector operator+(const Vector & right) const 
  {
    Vector result(*this);
    for(unsigned int i=0; i<n; i++)
      result[i]+=right[i];
    return result;
  }

  /**
   * Vector subtraction.
   * @param right The right operand of the minus sign.
   */
  Vector operator-(const Vector & right) const 
  {
    Vector result(*this);
    for(unsigned int i=0; i<n; i++)
      result[i]-=right[i];
    return result;
  }

  /**
	 * Fill the vector from the given pointer to an array.
	 * @param fillData An array that has at least 'n' elements.
	 */
	void set(T *fillData) {memcpy(this->data,fillData,n*sizeof(T));}

	/**
	 * Return the dimensionality of this vector.
	 */
	unsigned int size() {return n;}

  /**
   * Compute the l2 norm of the vector.
   */
  T l2Norm() {
    T result = 0;
    for(unsigned int i=0; i<n; i++)
      result+= data[i]*data[i];
    return sqrt(result);
  }

  /**
   * Make the vector unit sized (divide by the l2Norm).
   */
  void normalize() {
    T norm = this->l2Norm();
    for(unsigned int i=0; i<n; i++)
      data[i]/=norm;    
  }

private:

	/**
	 * Return this vector as a vnl_vector_ref object.
	 */
	vnl_vector_ref<T> getVnlVector() const {return vnl_vector_ref<T>(n, const_cast<double *>(this->data));}

  T data[n];
};


template<class T, unsigned int n>
inline Vector<T, n> operator*(const T & s, const Vector<T,n> & v) {return v*s;}

} //namespace lsqrRecipes

#endif //_VECTOR_H_
