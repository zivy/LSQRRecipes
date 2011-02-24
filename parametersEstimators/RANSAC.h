#ifndef _RANSAC_H_
#define _RANSAC_H_

#include <set>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits>

#include "ParametersEstimator.h"

/**
 * This class implements the RAndom SAmple Consensus (RANSAC) algorithm,
 * an algorithm for robust parameter estimation.
 * Given data containing outliers we estimate the model parameters using 
 * sub-sets of the original data:
 * 1. Choose the minimal subset from the data for computing the exact model 
 *    parameters.
 * 2. See how much of the input data agrees with the computed parameters.
 * 3. Goto step 1. This can be done up to (N choose m) times, where m is the 
 *    number of data objects required for an exact estimate and N is the total 
 *    number of data objects.
 * 4. Take the largest subset of objects which agreed on the parameters and 
 *    compute a least squares fit using them.
 * 
 * This is based on:
 * Fischler M.A., Bolles R.C., 
 * "Random Sample Consensus: A Paradigm for Model Fitting with Applications to 
 * Image Analysis and Automated Cartography", Communications of the ACM, 
 * Vol. 24(6), 1981.
 *
 * Hartely R., Zisserman A., "Multiple View Geometry in Computer Vision", 2001.
 *
 * The class template parameters are T - objects used for the parameter estimation 
 *                                      (e.g. Point2D in line estimation, 
 *                                            std::pair<Point2D,Point2D> in 
 *                                            homography estimation).
 *                                   S - type of parameter (e.g. double).                          
 *
 * @author: Ziv Yaniv (zivy@isis.georgetown.edu)
 *
 */

namespace lsqrRecipes {

template<class T, class S>
class RANSAC {

public:
	/**
	 * Estimate the model parameters using the RANSAC framework.
	 * @param parameters A vector which will contain the estimated parameters.
	 *                   If there is an error in the input then this vector will 
   *                   be empty.
	 *                   Errors are: 1. Less data objects than required for an 
   *                                  exact fit.
	 *                               2. The given data is in a singular 
   *                                  configuration (e.g. trying to fit a circle
	 *                                  to a set of colinear points).
   *                               3. The given parameter 
   *                                  desiredProbabilityForNoOutliers is not in 
   *                                 (0,1).
	 * @param paramEstimator An object which can estimate the desired parameters 
   *                       using either an exact fit or a least squares fit.
	 * @param data The input from which the parameters will be estimated.
	 * @param desiredProbabilityForNoOutliers The probability that at least one of
   *                                        the selected subsets doesn't contain 
   *                                        an outlier, must be in (0,1).
   * @param consensusSet Optional parameter. The consensus set, same length as 
   *                     the data vector, true if a data element belongs to the
   *                     set, false otherwise.
	 * @return Returns the percentage of data used in the least squares estimate.
	 */
	 static double compute(std::vector<S> &parameters, 
		                     ParametersEstimator<T,S> *paramEstimator , 
												 std::vector<T> &data, 
		                     double desiredProbabilityForNoOutliers,
                         std::vector<bool> *consensusSet=NULL);


	/**
	 * Estimate the model parameters using the maximal consensus set by going over
   * ALL possible subsets (brute force approach).
	 * Given: n -  data.size()
	 *        k - numForEstimate
	 * We go over all n choose k subsets       n!
	 *                                     ------------
	 *                                      (n-k)! * k!
	 * @param parameters A vector which will contain the estimated parameters.
	 *                   If there is an error in the input then this vector will 
   *                   be empty.
	 *                   Errors are: 1. Less data objects than required for an 
   *                                  exact fit.
	 *                               2. The given data is in a singular 
   *                                  configuration (e.g. trying to fit a circle
	 *                                  to a set of colinear points).
	 * @param paramEstimator An object which can estimate the desired parameters 
   *                       using either an exact fit or a least squares fit.
	 * @param data The input from which the parameters will be estimated.
	 * @param numForEstimate The number of data objects required for an exact fit.
   * @param consensusSet Optional parameter. The consensus set, same length as 
   *                     the data vector, true if a data element belongs to the
   *                     set, false otherwise.
	 * @return Returns the percentage of data used in the least squares estimate.
   *
	 * NOTE: This method should be used only when n choose k is small 
   *       (i.e. k or (n-k) are approximatly equal to n)
	 *
	 */
	 static double compute(std::vector<S> &parameters, 
		                     ParametersEstimator<T,S> *paramEstimator , 
                         std::vector<T> &data, std::vector<bool> *consensusSet=NULL);
private:

	 /**
	  * Compute n choose m  [ n!/(m!*(n-m)!)]. 
    * If choose(n,m)>std::numeric_limits<unsigned int>::max(), or there is an
    * overflow during the computations then we return 
    * std::numeric_limits<unsigned int>::max(), otherwise the correct value
    * is returned.
		*/
	static unsigned int choose(unsigned int n, unsigned int m);

	static void computeAllChoices(ParametersEstimator<T,S> *paramEstimator, 
                                std::vector<T> &data,
												        bool *bestVotes, bool *curVotes, 
                                int &numVotesForBest, int startIndex, int k, 
                                int arrIndex, int *arr);

	static void estimate(ParametersEstimator<T,S> *paramEstimator, 
                       std::vector<T> &data, bool *bestVotes, bool *curVotes, 
                       int &numVotesForBest, int *arr);

	 class SubSetIndexComparator {
		 private:
			int length;
		 public:
			SubSetIndexComparator(int arrayLength) : length(arrayLength){}
			bool operator()(const int *arr1, const int *arr2) const {
        for(int i=0; i<this->length; i++) {
					if(arr1[i] < arr2[i])
						return true;
          else if(arr1[i] > arr2[i]) 
            return false;
        }
				return false;			
			}
		};
};

} //namespace lsqrRecipes

#include "RANSAC.txx"

#endif //_RANSAC_H_
