#ifndef _DENSE_LINEAR_EQUATION_SYSTEM_PARAMETERS_ESTIMATOR_H_ 
#define _DENSE_LINEAR_EQUATION_SYSTEM_PARAMETERS_ESTIMATOR_H_

#include<cstring>
#include "copyright.h"

#include "ParametersEstimator.h"


namespace lsqrRecipes {


/**
 * Template defining the n+1 augmented row structure of the matrix equation 
 * system Ax=b, [a_0,a_1,...a_{n-1},b]. Valid types for row entries are the 
 * numeric types that represent real numbers (float, double).
 *
 * @author: Ziv Yaniv
 */
template <class T, unsigned int n>
class AugmentedRow {
public:
            //enable compile time access to the row dimension, ugly but necessary
  enum {dimension=n};

  /**
	 * Output the augmented row data to the given output stream as the following string:
	 *[a_0,a_1,...a_{n-1},b]
	 */
  friend std::ostream &operator<<(std::ostream& output, const AugmentedRow &r) {
		  output<<"[ ";
		  for(unsigned int i=0; i<n; i++) 
		    output<<r.aValues[i]<<", ";
      output<<r.bValue<<" ]";
		  return output;
	  }

   /**
    * Assignment operator.
    */
  AugmentedRow<T,n> &operator =(const AugmentedRow<T,n> &other) {
    memcpy(this->aValues, other.aValues, n*sizeof(T)); this->bValue = other.bValue;
    return *this;
  }

	 /**
	  * Access to the entries of the augmented row. No bounds checking is performed.
    * Valid indexes are in [0,n], corresponding to [a_0...a_{n-1},b].
	  */
  T & operator[](unsigned int index) {return index==n ? this->bValue : this->aValues[index];}
  const T & operator[](unsigned int index) const {return index==n ? this->bValue : this->aValues[index];}

   /**
	  * Default constructor, sets the augmented row to zero.
	  */
  AugmentedRow() {
    memset(this->aValues,0,n*sizeof(T));
    this->bValue = static_cast<T>(0.0);
  }

   /**
	  * Construct an augmented row from the given pointer to an array.
	  * @param fillData An array that has at least 'n+1' elements.
	  */
  AugmentedRow(T *fillData) {
    memcpy(this->aValues,fillData,n*sizeof(T));
    this->bValue = fillData[n];
  }

   /**
	  * Construct an augmented row from the given pointer to an array and b value.
	  * @param fillData An array that has at least 'n' elements.
    * @param bData Set the augmented element to this value.
	  */
  AugmentedRow(T *fillData, T bData) {
    memcpy(this->aValues,fillData,n*sizeof(T));
    this->bValue = bData;
  }

   /**
	  * Copy constructor.
	  * @param other The augmented row we copy.
	  */
  AugmentedRow(const AugmentedRow<T,n> &other) {
    memcpy(this->aValues,other.aValues,n*sizeof(T));
    this->bValue = other.bValue;
  }
	
   /**
	  * Fill the augmented row from the given pointer to an array.
	  * @param fillData An array that has at least 'n+1' elements.
	  */
  void set(T *fillData) {
    memcpy(this->aValues,fillData,n*sizeof(T)); 
    this->bValue = fillData[n];
  }

   /**
	  * Fill the augmented row from the given pointer to an array and b value.
	  * @param fillData An array that has at least 'n' elements.
    * @param bData Set the augmented element to this value.
	  */
  void set(T *fillData, T bData) {
    memcpy(this->aValues,fillData,n*sizeof(T));
    this->bValue = bData;
  }

   /**
	  * Get the augmented row data.
	  * @param aValues An array that has at least 'n' elements.
    * @param bValue A scalar of the same type as the array.
	  */
  void get(T *aValues, T &bValue) {
    memcpy(aValues, this->aValues,n*sizeof(T)); 
    bValue = this->bValue;
  }

   /**
	  * Get the augmented row data.
	  * @param aValues An array that has at least 'n+1' elements.
	  */
  void get(T *aValues) {
    memcpy(aValues, this->aValues,n*sizeof(T));
    aValues[n] = this->bValue;
  }

  /**
	  * Return the dimensionality of the augmented row.
	  */
  unsigned int size() {return n+1;}
 
private:
  T aValues[n];
  T bValue;
};


/**
 * This class estimates the parameters/solution of a dense overdetermined linear 
 * equation system in a least squares sense, x = argmin_x \|Ax-b\|. 
 * The solution is given by the normal equations: x = (A^TA)^{-1}A^Tb
 * We provide the unique solution if it exists (rank(A)=colNum), otherwise we
 * return an empty vector.
 * 
 * @author: Ziv Yaniv
 *
 */
template<class T, unsigned int n>  //n is the number of columns in A
class DenseLinearEquationSystemParametersEstimator : public ParametersEstimator< AugmentedRow<T, n>, T> {
public:



  /**
	 * Object constructor.
	 * @param delta An equation a^Tx= b is consistent with the solution x if
   *              abs(a^Tx - b)< delta .
   */
	DenseLinearEquationSystemParametersEstimator(T delta);

	/**
	 * Compute the solution to the equation system defined by the given augmented 
   * rows.
	 * @param data A vector containing n augmented rows of the form [a_0,a_1,...a_{n-1},b].
	 * @param parameters This vector is cleared and then filled with the computed 
   *                   parameters, The solution to the equation system:
   *                   [a_{0,0} a_{0,1} ... a_{0,n-1}] [x_0]      [b_0]
   *                   [  :      :           :       ] [ : ]    = [ : ] 
   *                   [a_{n,0} a_{n,1} ... a_{n,n-1}] [x_{n-1}]  [b_{n-1}]
   *                 
   *                   If the data contains less than n augmented rows, or 
   *                   A is not invertible then the resulting parameters 
   *                   vector is empty (size = 0).
	 */
  virtual void estimate(std::vector< AugmentedRow<T, n> *> &data, 
                        std::vector<T> &parameters);
  virtual void estimate(std::vector< AugmentedRow<T, n> > &data, 
                        std::vector<T> &parameters);

	/**
	 * Compute a least squares solution to the equation system defined by the given
   * augmented rows.
	 *
	 * @param data The solution minimizes the least squares error of the equation 
   *             system.
	 * @param parameters This vector is cleared and then filled with the computed 
   *                   parameters, The solution to the over determined equation 
   *                   system: Ax = b
   *                 
   *                   If the data contains less than n augmented rows, or 
   *                   rank(A)!= n (no unique solution) then the resulting 
   *                   parameters vector is empty (size = 0).	 
   */
  virtual void leastSquaresEstimate(std::vector< AugmentedRow<T, n> *> &data, 
                                    std::vector<T> &parameters);
  virtual void leastSquaresEstimate(std::vector< AugmentedRow<T, n> > &data, 
                                    std::vector<T> &parameters);

	/**
	 * Return true if the error for the specific equation is smaller than 'delta' 
   * (see constructor).
	 * @param parameters The potential solution [x_0,...,x_{n-1}].
	 * @param data Check that the error, |dot(a,x) - b|, is less than 'delta'.
   *
   * Note:If the parameters vector is too short the method will throw an 
   *      exception as it tries to access a non existing entry. 
   *
	 */
	virtual bool agree(std::vector<T> &parameters, AugmentedRow<T, n> &data);

	/**
	 * Change the tolerance defining whether an equation is consistent with a 
   * given solution.
	 * @param delta An equation a^Tx= b is consistent with the solution x if
   *              abs(a^Tx - b)< delta.
	 */
	void setDelta(T delta) {this->delta = delta;}

  /**
   * Utility method for transforming an equation system represented by a vnl
   * matrix and vector into the representation used by our class, augmented rows.
   */
  static void getAugmentedRows(vnl_matrix<T> &A, vnl_vector<T> &b, 
                               std::vector< AugmentedRow<T, n> > &rows);

private:
    //given solution x and equation a^Tx=b, if abs(a^Tx - b) < delta then the 
    //equation is consistent with the solution
	double delta; 
};

} //namespace lsqrRecipes

#include "DenseLinearEquationSystemParametersEstimator.hxx" //the implementation is in this file

#endif //_DENSE_LINEAR_EQUATION_SYSTEM_PARAMETERS_ESTIMATOR_H_
