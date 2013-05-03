namespace lsqrRecipes {

template<class T, class S>
double RANSAC<T,S>::compute(std::vector<S> &parameters, 
                            ParametersEstimator<T,S> *paramEstimator , 
                            std::vector<T> &data,  		                     
                            double desiredProbabilityForNoOutliers,
                            std::vector<bool> *consensusSet)
{
  unsigned int numDataObjects = static_cast<unsigned int>(data.size());
  unsigned int numForEstimate = paramEstimator->numForEstimate();

		      //there are less data objects than the minimum required for an 
		      //exact fit, or
	        //desiredProbabilityForNoOutliers is not in (0.0,1.0)
  if(numDataObjects < numForEstimate || 
     desiredProbabilityForNoOutliers>=1.0 || 
     desiredProbabilityForNoOutliers<=0.0) 
    return 0;

  std::vector<T *> exactEstimateData;
  std::vector<T *> leastSquaresEstimateData;
  std::vector<S> exactEstimateParameters;
  unsigned int i, k, l, m, numVotesForBest, numVotesForCur, maxIndex;
  int j;
  unsigned int numTries;
                      //true if data[i] agrees with the best model, otherwise false
  bool *bestVotes = new bool[numDataObjects]; 
                 //true if data[i] agrees with the current model, otherwise false
  bool *curVotes = new bool[numDataObjects];
               //true if data[i] is NOT chosen for computing the exact fit, otherwise false
  bool *notChosen = new bool[numDataObjects]; 
  SubSetIndexComparator subSetIndexComparator(numForEstimate);
  std::set<int *, SubSetIndexComparator > chosenSubSets(subSetIndexComparator);
  int *curSubSetIndexes;
  double numerator = log(1.0-desiredProbabilityForNoOutliers);
  double denominator;
	              //allTries is either the correct value or if there is an overflow
	              //during the computation it is set to the maximal value
	              //for unsigned int
  unsigned int allTries = choose(numDataObjects,numForEstimate);

  parameters.clear();
  srand((unsigned)time(NULL)); //seed random number generator

  numVotesForBest = 0; //initalize with 0 so that the first computation which gives any type of fit will be set to best
  numTries = allTries; //initialize with the number of all possible subsets
	    
  for(i=0; i<numTries; i++) {
				  //randomly select data for exact model fit ('numForEstimate' objects).
    std::fill(notChosen,notChosen+numDataObjects, true); 
    curSubSetIndexes = new int[numForEstimate];
		
    exactEstimateData.clear();
		
    maxIndex = numDataObjects-1; 
    for(l=0; l<numForEstimate; l++) {
		                             //selectedIndex is in [0,maxIndex]
      int selectedIndex = (int)(((float)rand()/(float)RAND_MAX)*maxIndex + 0.5);
      for(j=-1,k=0; k<numDataObjects && j<selectedIndex; k++) {
        if(notChosen[k])
          j++;
      }
			k--;
			exactEstimateData.push_back(&(data[k]));
			notChosen[k] = false;
			maxIndex--;
		}
                 //get the indexes of the chosen objects so we can check that this sub-set hasn't been
		             //chosen already
		for(l=0, m=0; m<numDataObjects; m++) {
			if(!notChosen[m]) {
				curSubSetIndexes[l] = m+1;
				l++;
			}
		}

                 //check that the sub-set just chosen is unique
		std::pair< typename std::set<int *, SubSetIndexComparator >::iterator, bool> res = chosenSubSets.insert(curSubSetIndexes);
		 
		if(res.second == true) { //first time we chose this sub set
		
				                 //use the selected data for an exact model parameter fit
      paramEstimator->estimate(exactEstimateData,exactEstimateParameters);
			                  //selected data is a singular configuration (e.g. three colinear points for 
			                  //a circle fit)
			if(exactEstimateParameters.size() == 0)
				continue;
					     //see how many agree on this estimate
			numVotesForCur = 0;
      std::fill(curVotes,curVotes+numDataObjects, false);			
                                //continue checking data until there is no chance of getting a larger consensus set 
                                //or all the data has been checked              
			for(m=0; m<numDataObjects && numVotesForBest-numVotesForCur<numDataObjects-m+1; m++) {
				if(paramEstimator->agree(exactEstimateParameters, data[m])) {
					curVotes[m] = true;
					numVotesForCur++;
				}
			}                            //found a larger consensus set?
			if(numVotesForCur > numVotesForBest) {
				numVotesForBest = numVotesForCur;				
        std::copy(curVotes, curVotes+numDataObjects, bestVotes);
                      //all data objects are inliers, terminate the search
			  if(numVotesForBest == numDataObjects)
          i=numTries;                
        else {  //update the estimate of outliers and the number of iterations we need				           		  				
				  denominator = log(1.0- pow((double)numVotesForCur/(double)numDataObjects, (double)(numForEstimate)));
				  numTries = (int)(numerator/denominator + 0.5);
					            //there are cases when the probabilistic number of tries is greater than all possible sub-sets
				  numTries = numTries<allTries ? numTries : allTries;
        }
      }
		}
		else {  //this sub set already appeared, release memory
			delete [] curSubSetIndexes;			
		}
	}
                 
	     //release the memory
	typename std::set<int *, SubSetIndexComparator >::iterator it = chosenSubSets.begin();
	typename std::set<int *, SubSetIndexComparator >::iterator chosenSubSetsEnd = chosenSubSets.end();
	while(it!=chosenSubSetsEnd) {
		delete [] (*it);
		it++;
	}
	chosenSubSets.clear();

	             //compute the least squares estimate using the largest sub set
	if(numVotesForBest > 0) {
		for(m=0; m<numDataObjects; m++) {
			if(bestVotes[m])
				leastSquaresEstimateData.push_back(&(data[m]));
		}
		if(consensusSet!=NULL) {
		  consensusSet->clear();
		  consensusSet->insert(consensusSet->begin(), bestVotes, bestVotes+numDataObjects);
		}
		paramEstimator->leastSquaresEstimate(leastSquaresEstimateData,parameters);
	}
	delete [] bestVotes;
	delete [] curVotes;
	delete [] notChosen;

	return (double)numVotesForBest/(double)numDataObjects;
}

/*****************************************************************************/

template<class T, class S>
double RANSAC<T,S>::compute(std::vector<S> &parameters, 
													  ParametersEstimator<T,S> *paramEstimator , 
												    std::vector<T> &data,
												    std::vector<bool> *consensusSet)
{
  unsigned int numForEstimate = paramEstimator->numForEstimate();
	std::vector<T *> leastSquaresEstimateData;
	unsigned int numDataObjects = static_cast<unsigned int>(data.size());
	unsigned int numVotesForBest = 0;
	int *arr = new int[numForEstimate];
            //true if data[i] agrees with the current model, otherwise false
	bool *curVotes = new bool[numDataObjects];  
              //true if data[i] agrees with the largest consensus model
	bool *bestVotes = new bool[numDataObjects];  
	
	parameters.clear();

		      //there are less data objects than the minimum required for an exact fit
	if(numDataObjects < numForEstimate) 
		return 0;

	computeAllChoices(paramEstimator,data,
										bestVotes, curVotes, numVotesForBest, 0, numForEstimate, 0, arr);
	
	   //compute the least squares estimate using the largest sub set
	if(numVotesForBest > 0) {
		for(unsigned int j=0; j<numDataObjects; j++) {
			if(bestVotes[j])
				leastSquaresEstimateData.push_back(&(data[j]));
		}
		if(consensusSet!=NULL) {
		  consensusSet->clear();
		  consensusSet->insert(consensusSet->begin(), bestVotes, bestVotes+numDataObjects);
		}
		paramEstimator->leastSquaresEstimate(leastSquaresEstimateData,parameters);
	}

	delete [] arr;
	delete [] bestVotes;
	delete [] curVotes;	

	return (double)numVotesForBest/(double)numDataObjects;
}

/*****************************************************************************/

template<class T, class S>
void RANSAC<T,S>::computeAllChoices(ParametersEstimator<T,S> *paramEstimator, std::vector<T> &data,
																		bool *bestVotes, bool *curVotes, unsigned int &numVotesForBest, int startIndex, int k, int arrIndex, int *arr)
{
	              //we have a new choice of indexes
  if(k==0) {
		estimate(paramEstimator, data, bestVotes, curVotes, numVotesForBest, arr);
    return;
  }
	       //continue to recursivly generate the choice of indexes
  int endIndex = data.size()-k;
  for(int i=startIndex; i<=endIndex; i++) {
    arr[arrIndex] = i;
    computeAllChoices(paramEstimator, data, bestVotes, curVotes, numVotesForBest,
			                i+1, k-1, arrIndex+1, arr);
  }

}
/*****************************************************************************/

template<class T, class S>
void RANSAC<T,S>::estimate(ParametersEstimator<T,S> *paramEstimator, std::vector<T> &data, 
											     bool *bestVotes, bool *curVotes, unsigned int &numVotesForBest, int *arr)
{
	std::vector<T *> exactEstimateData;
	std::vector<S> exactEstimateParameters;
	unsigned int numDataObjects;
	unsigned int numVotesForCur;
	unsigned int j;

	numDataObjects = static_cast<unsigned int>(data.size());
  std::fill(curVotes,curVotes+numDataObjects, false);
	numVotesForCur=0;

  unsigned int numForEstimate = paramEstimator->numForEstimate();

	for(j=0; j<numForEstimate; j++)
		exactEstimateData.push_back(&(data[arr[j]]));
	paramEstimator->estimate(exactEstimateData,exactEstimateParameters);
	                     //singular data configuration
	if(exactEstimateParameters.size()==0)
		return;

	for(j=0; j<numDataObjects; j++) {
		if(paramEstimator->agree(exactEstimateParameters, data[j])) {
			curVotes[j] = true;
			numVotesForCur++;
		}
	}
	if(numVotesForCur > numVotesForBest) {
		numVotesForBest = numVotesForCur;
    std::copy(curVotes, curVotes+numDataObjects, bestVotes);		
	}
}

/*****************************************************************************/

template<class T, class S>
unsigned int RANSAC<T,S>::choose(unsigned int n, unsigned int m)
{
	double denominatorEnd, numeratorStart, numerator,denominator, i, result; 
        //perform smallest number of multiplications
	if((n-m) > m) {
		numeratorStart = n-m+1;
		denominatorEnd = m;
	}
	else {
		numeratorStart = m+1;
		denominatorEnd = n-m;
	}
	
	for(i=numeratorStart, numerator=1; i<=n; i++)
		numerator*=i;
	for(i=1, denominator=1; i<=denominatorEnd; i++)
		denominator*=i;	
	result = numerator/denominator;
	
	         //check for overflow both in computation and in result	         
	if(denominator>std::numeric_limits<double>::max() || 
	   numerator>std::numeric_limits<double>::max() || 
	   static_cast<double>(std::numeric_limits<unsigned int>::max())<result )	
	  return std::numeric_limits<unsigned int>::max();
	else 
	  return static_cast<unsigned int>(result);   
}

} //namespace lsqrRecipes
