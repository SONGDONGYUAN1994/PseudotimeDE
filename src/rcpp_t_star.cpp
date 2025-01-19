
/***
 * A collection of functions that implement the main functionality of the improved
 * version of algorithms based in
 *
 * Weihs, Luca, Mathias Drton, and Dennis Leung. "Efficient Computation of the
 * Bergsma-Dassios Sign Covariance." arXiv preprint arXiv:1504.00964 (2015).
 */

#include<RcppArmadilloExtensions/sample.h>
// #include"red_black_tree.h"
#include<stdio.h>
#include<ctype.h>
#include<cmath>
#include<algorithm>
#include<iostream>
#include <cstdint>
#include <vector>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

// The following functions are used as input to initalize a red black tree whose
// nodes contain doubles. Knowing why these functions exist is likely
// unimportant to you, if otherwise then please read the red_black_tree.cpp file
// to see where they are used.
void DoubDest(void* a) { free((double*)a); }
int DoubComp(const void* a, const void* b) {
	if (*(double*)a > *(double*)b) {
    return(1);
	}
	if (*(double*)a < *(double*)b) {
    return(-1);
	}
	return(0);
}
void DoubPrint(const void* a) { Rprintf("%f", *(double*)a); }
void InfoPrint(void* a) { }
void InfoDest(void *a) { }

/***
 * A simple bubble sort implementation, this is only ever used on a vec of size
 * 4 in the following functions and hence its inefficiency is unimportant.
 * @param vec the array to sort
 * @param n the number of elements in the vec array
 */
void bubbleSort(long double * vec, int n) {
	int i,j;
	long double tmp;
	for(i = n - 1; i > 0; i--) {
		for(j = 0; j < i; j++) {
			if(vec[j + 1] < vec[j]) {
				tmp = vec[j + 1];
				vec[j + 1] = vec[j];
				vec[j] = tmp;
			}
		}
	}
}

/***
 * For an input vector x = (x_1,...,x_n) returns a new vector r of length n
 * with r[i] equal to the rank of x_i in x. So if x == (4,2,2,3,2,3) then
 * r == (3, 1, 1, 2, 1, 2).
 * @param x an arma::vec to get the ranks of
 * @return the vector of ranks
 */
arma::uvec vecToRanks(const arma::vec& x) {
  if (x.n_elem == 0) {
    return(arma::zeros<arma::uvec>(0));
  }
  arma::uvec xSortedInds = arma::sort_index(x);
  arma::uvec ranks = arma::zeros<arma::uvec>(x.n_elem);
  arma::uvec repeats = arma::zeros<arma::uvec>(x.n_elem);
  int k = 1;
  double last = x[xSortedInds[0]];
  for (unsigned int i = 0; i < x.n_elem; i++) {
    double current = x[xSortedInds[i]];
    if (current != last) {
      k++;
      last = current;
    }
    repeats[i] = k;
  }
  for (unsigned int i = 0; i < xSortedInds.n_elem; i++) {
    ranks[xSortedInds[i]] = repeats[i];
  }
  return ranks;
}

/***
 * For two input vectors of ranks xRanks, yRanks (i.e. taking values in the
 * positive integers and including all integers between 1 and their maximum
 * value) returns the matrix A where A[i,j] is the number times and element of
 * xRanks, and the corresponding element of yRanks, are less than or equal to i
 * and j respectively. Hence (xRanks[s],yRanks[s]) contributes to A[i,j] if and
 * only if x[s] <= i and y[s] <= i.
 * @param xRanks a vector of ranks.
 * @param yRanks a vector of ranks (must be the same length as x).
 * @return the "less than or equal" to matrix.
 */
arma::Mat<long long> ranksToLeqMat(const arma::uvec& xRanks, const arma::uvec& yRanks) {
  int xMax = xRanks.max();
  int yMax = yRanks.max();
  arma::Mat<long long> leqMat = arma::zeros<arma::Mat<long long>>(xMax + 1, yMax + 1);
  for (unsigned int i = 0; i < xRanks.n_elem; i++) {
    leqMat(xRanks[i], yRanks[i]) += 1;
  }
  for (int i = 1; i <= xMax; i++) {
    for (int j = 1; j <= yMax; j++) {
      leqMat(i,j) =
        leqMat(i,j - 1) + leqMat(i - 1,j) - leqMat(i - 1, j - 1) + leqMat(i,j);
    }
  }
  return leqMat;
}

/***
 * Computes a version of the matrix defined at the top of page 5 of Heller and
 * Heller (2016) (arxiv:1605.08732). In particular, we let the rank by always
 * <= rather than < in the x case. This makes the notation nice in
 * TStarHellerAndHellerRCPP.
 * @param leqMat a leqMat as created by ranksToLeqMat.
 * @return the unique count matrix.
 */

arma::Mat<long long> leqMatToUniqueCountMat(const arma::Mat<long long>& leqMat) {
  arma::Mat<long long> uCountMat = arma::zeros<arma::Mat<long long>>(leqMat.n_rows, leqMat.n_cols);
  for (unsigned int i = 1; i < leqMat.n_rows; i++) {
    for (unsigned int j = 1; j < leqMat.n_cols; j++) {
      int numEqYAndLessOrEqX = leqMat(i, j) - leqMat(i, j - 1);
      uCountMat(i,j) = uCountMat(i, j - 1) +
        (numEqYAndLessOrEqX * (numEqYAndLessOrEqX - 1)) / 2;
    }
  }
  return uCountMat;
}

/***
 * Given an arma::uvec x and a set of indices inds, creates a new vector whose
 * entries are those of x indexed by inds. So if
 * x = (1,2,3)
 * inds = (0,0,2)
 * then we would return the vector (1,1,3).
 * @param x vector be be indexed into
 * @param inds the indices
 * @return a new vector.
 */
arma::uvec indexUvec(const arma::uvec& x, const arma::uvec& inds) {
  arma::uvec newVec = arma::zeros<arma::uvec>(inds.n_elem);
  for (unsigned int i = 0; i < newVec.n_elem; i++) {
    newVec[i] = x[inds[i]];
  }
  return newVec;
}

/***
 * Function implementing an improved computation of the t* U-statistic
 * based in
 *
 * Heller, Yair and Heller, Ruth. "Computing the Bergsma Dassios sign-
 * covariance." arXiv preprint arXiv:1605.08732 (2016).
 *
 * @param x a arma::vec of values
 * @param y a arma::vec of values of the same length as x
 * @return the U-statistic computed upon the two input vectors.
 */
// [[Rcpp::export]]

double TStarHellerAndHellerRCPP(const arma::vec& x, const arma::vec& y) {
  arma::uvec xRanks = vecToRanks(x);
  arma::uvec yRanks = vecToRanks(y);
  arma::Mat<long long> leqMat = ranksToLeqMat(xRanks, yRanks);
  arma::Mat<long long> uCountMat = leqMatToUniqueCountMat(leqMat);
  arma::uvec xOrder = arma::sort_index(xRanks);
  xRanks = indexUvec(xRanks, xOrder);
  yRanks = indexUvec(yRanks, xOrder);

	
  ///////////////////////////////////////////////////
  ///////////////////////////////////////////////////
  int ymax = yRanks.max();
  int xmax = xRanks[xRanks.n_elem - 1];

  long long numCon = 0;
  long long numDis = 0;
  std::vector<long long> yRankCount(ymax + 1, 0);
  for (int i = 1; i <= ymax; i++) {
    yRankCount[i] = leqMat(xmax, i) - leqMat(xmax, i - 1);
  }
  
  for (unsigned int i = 0; i < xRanks.n_elem - 1; i++) {
    yRankCount[yRanks[i]] -= 1;
    int xRankMin = xRanks[i];
    for (int cnt = 1; cnt <= ymax; cnt++){ // loop over the unique values of y
      int yRankMin = std::min(int(yRanks[i]), cnt);
      int yRankMax = std::max(int(yRanks[i]), cnt);

      long long bot = leqMat(xRankMin - 1, yRankMin - 1);
      long long mid = (yRankMin == yRankMax) ? 0 :
        leqMat(xRankMin - 1, yRankMax - 1) - leqMat(xRankMin - 1, yRankMin);
      long long top = leqMat(xRankMin - 1, ymax) -
                 leqMat(xRankMin - 1, yRankMax);
      long long eqMin = leqMat(xRankMin - 1, yRankMin) -
        leqMat(xRankMin - 1, yRankMin - 1);
      long long eqMax = leqMat(xRankMin - 1, yRankMax) -
        leqMat(xRankMin - 1, yRankMax - 1);

      numCon += yRankCount[cnt] * (top * (top - 1) / 2.0 + bot * (bot - 1) / 2.0);

      if (yRankMin != yRankMax) {
        numDis += yRankCount[cnt] * (top * (mid + eqMin + bot) + bot * (mid + eqMax) +
          eqMin * (mid + eqMax) +
          eqMax * mid + mid * (mid - 1) / 2.0);
        numDis -= yRankCount[cnt] * (uCountMat(xRankMin - 1, yRankMax - 1) -
          uCountMat(xRankMin - 1, yRankMin));
      }
      
    }
  }

  int n = xRanks.n_elem;
  long long c = 16 * numCon - 8 * numDis;
  double d = (c < 0) ? -1 : 1;
  return d * expl(logl(d * c) - (logl(n) + logl(n - 1) +
                               logl(n - 2) + logl(n - 3)));

}
