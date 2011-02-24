#ifndef _EPSILON_H_
#define _EPSILON_H_

#include "copyright.h"

/**
 * How close is close enough? 
 * File containing all the relevant choices for epsilon.
 *
 * @author: Ziv Yaniv (zivy@isis.georgetown.edu)
 */

namespace lsqrRecipes {

                //I chose to set it to be DBL_EPSILON (see float.h) which is the smallest positive number x, 
               //such that x + 1.0 is not equal to 1.0, another option is FLT_EPSILON==1.192092896e–07F which is
              //the same for single precision numbers               
#define EPS 2.220446049250313e-016

} //namespace lsqrRecipes

#endif //_EPSILON_H_
