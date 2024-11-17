/***
 * Copyright (C) 2015 Luca Weihs
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/***
 * A collection of functions that allow for the computation of the asymptotic
 * distributions described in Nandy, Weihs, Drton (2016)
 * <http://arxiv.org/abs/1602.04387>.
 */

#include<RcppArmadillo.h>
#include<algorithm>
#include<cmath>
#include<queue>
#include "AsymMixedCdfIntegrandEvaluator.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

using namespace Rcpp;

/***
 * Used in the numericalCfInversion function. Doubles the range of values
 * considered by the numerical integration; that is, if we were previously
 * integrating in the range [-T, T] this extends the range to [-2*T, 2*T]).
 * If we using n points in the riemann sum prior to running doubleWidth
 * then we use 2*n points afterwards.
 */
void doubleWidth(arma::vec& positions, arma::vec& values,
                 IntegrandEvaluator& intEval, double x,
                 double integrandError) {
  int oldSize = values.size();
  double oldMaxPosition = positions(positions.size() - 1);
  double newMaxPosition = 2 * oldMaxPosition;
  double intervalSize = (newMaxPosition - oldMaxPosition) / oldSize;
  
  positions.resize(oldSize * 2);
  values.resize(oldSize * 2);
  for (int i = oldSize; i < 2 * oldSize; i++) {
    positions(i) = positions(i - 1) + intervalSize;
    values(i) = 2.0 * intEval.integrand(x, positions(i), integrandError).real();
  }
}

/***
 * Used in the numericalCfInversion function. Bisects every interval over which
 * the riemann integration occurs; that is, this makes the grid of points we
 * evaluate over finer.
 */
void bisect(arma::vec& positions, arma::vec& values,
            IntegrandEvaluator& intEval, double x,
            double integrandError) {
  int oldSize = values.size();
  
  positions.resize(oldSize * 2 - 1);
  values.resize(oldSize * 2 - 1);
  for (int i = oldSize - 1; i > 0; i--) {
    positions(2 * i) = positions(i);
    values(2 * i) = values(i);
  }
  for (unsigned int i = 1; i < positions.size(); i += 2) {
    positions(i) = (positions(i + 1) + positions(i - 1)) / 2.0;
    values(i) = 2.0 * intEval.integrand(x, positions(i), integrandError).real();
  }
}

/***
 * Used in the numericalCfInversion function. Takes a vector of ordered
 * positions (points) in R along with some collection of values associated to
 * these points and then computes a riemann sum of the values using the points
 * as a partition. Here the size of an interval associated with a value is
 * determined by the distance between the point's position and the next point's
 * position.
 */
double riemannIntegrate(const arma::vec& positions, const arma::vec& values) {
  if (positions(0) != 0 && positions.size() >= 2) {
    stop("riemannIntegrate expects the first position to be 0 and"
           " there must be at least 2 positions.");
  }
  double sum = 0;
  
  // Special case handling of 0th position
  double bisectVal = positions(1) / 2;
  sum += values(0) * (2 * bisectVal);
  
  // Looping through position 1,...,n - 2
  for (unsigned int i = 1; i <= positions.size() - 2; i++) {
    double lastBisectVal = bisectVal;
    bisectVal = (positions(i + 1) + positions(i)) / 2;
    double intervalWidth = (bisectVal - lastBisectVal);
    sum += values(i) * intervalWidth;
  }
  
  // Special case handling of last position
  int lastInd = values.size() - 1;
  sum += values(lastInd) * (2 * (positions(lastInd) - bisectVal));
  
  return sum;
}

/***
 * Takes an integrand evaluator which corresponds to the integrand of a
 * numerical characteristic function inversion and then performs the inversion
 * using that integrand evaluator. The integration starts on the interval [-T,T]
 * but this interval grows to ensure that good convergence properties hold. We
 * assume that probability density or cumulative distribution function we are
 * trying to determine by inverting is supported on the positive real numbers
 * and thus this function will return 0 if x <= 0.
 *
 * @param intEval an IntegrandEvaluator object which provides the integrand used
 *        during the numerical integration.
 * @param x the value at which we wish to evaluate the inverted function.
 * @param T used to set the initial interval [-T,T] over which the integration
 *        occurs. Setting this intelligently can greatly improve performance.
 * @param convCrit a convergence criterion, the numerical integration will stop
 *        when within changes are less than convCrit.
 * @param maxIter the maximum number of iterations in the algorithm before
 *        stopping if the convCrit is not reached. If maxIter is reached then a
 *        warning is printed. If maxIter is < 5 then maxIter is set to 5.
 */
double numericalCfInversion(IntegrandEvaluator& intEval, double x, double T,
                            double convCrit, int maxIter) {
  // We implicitly assume that we work only with CDFs that are supported
  // strictly on the non-negative axis. As such we return 0 when x <= 0.
  if (x <= 0) {
    return 0;
  }
  double integrandError = 0.1; // convCrit * .0001
  
  int numInts = 5;
  double intWidth = T / numInts;
  arma::vec positions(numInts);
  arma::vec values(numInts);
  
  for (int j = 0; j < numInts; j++) {
    positions(j) = j * intWidth;
    if (j == 0) {
      values(j) = intEval.integrand(x, positions(j), integrandError).real();
    } else {
      values(j) = 2.0*intEval.integrand(x, positions(j), integrandError).real();
    }
  }
  
  double oldIntVal = riemannIntegrate(positions, values);
  
  bisect(positions, values, intEval, x, integrandError);
  double intVal = riemannIntegrate(positions, values);
  double bisectChange = std::fabs(static_cast<double>(oldIntVal - intVal)) + convCrit + 1;
  oldIntVal = intVal;
  
  doubleWidth(positions, values, intEval, x, integrandError);
  double widthChange = std::fabs(static_cast<double>(oldIntVal - intVal)) + convCrit + 1;
  
  int k = 0;
  while (k < 1 || (std::max(bisectChange, widthChange) >= convCrit && k < maxIter)) {
    oldIntVal = intVal;
    if (bisectChange > widthChange) {
      bisect(positions, values, intEval, x, integrandError);
      intVal = riemannIntegrate(positions, values);
      bisectChange = std::fabs(static_cast<double>(oldIntVal - intVal));
    } else {
      doubleWidth(positions, values, intEval, x, integrandError);
      intVal = riemannIntegrate(positions, values);
      widthChange = std::fabs(static_cast<double>(oldIntVal - intVal));
    }
    k++;
  }
  
  if (k == maxIter) {
    Rcpp::warning("Max iterations reached, cannot guarentee convergence.\n");
  }
  
  return intVal;
}

/***
 * A simple function used to bound values to be within [0,1].
 */
double boundInZeroOne(double x) {
  return std::min(std::max(x, 0.0), 1.0);
}

//' Eigenvalues for discrete asymptotic distribution
 //'
 //' Computes the eigenvalues needed to determine the asymptotic distributions
 //' in the mixed/discrete cases. See Nandy, Weihs, and Drton (2016)
 //' <http://arxiv.org/abs/1602.04387> for more details.
 //'
 //' @export
 //'
 //' @param p a vector of probabilities that sum to 1.
 //'
 //' @return the eigenvalues associated to the matrix generated by p
 // [[Rcpp::export]]
 arma::vec eigenForDiscreteProbs(arma::vec p) {
   arma::vec cdf(p.size());
   arma::vec q1(p.size());
   cdf[0] = p[0];
   for (unsigned int i = 1; i < p.size(); i++) {
     cdf[i] = cdf[i - 1] + p[i];
   }
   q1[0] = p[0] * (1 - cdf[0]);
   for (unsigned int i = 1; i < p.size(); i++) {
     q1[i] = q1[i-1] + p[i] * (1 - cdf[i]);
   }
   
   arma::mat symMat(p.size(), p.size());
   for(unsigned int i = 0; i < p.size(); i++) {
     for(unsigned int j = i; j < p.size(); j++) {
       if (i == 0) {
         symMat(i,j) = (1 - cdf[j]) * (1 - cdf[j]);
       } else {
         symMat(i,j) = (1 - cdf[j]) * (1 - cdf[j]) + cdf[i-1] * cdf[i-1];
       }
       
       if (i != j) {
         symMat(i,j) -= cdf[i] * (1 - cdf[i]) + (q1[j-1] - q1[i]);
         symMat(j,i) = symMat(i,j);
         symMat(j,i) *= std::sqrt(p[i] * p[j]);
       }
       symMat(i,j) *= std::sqrt(p[i] * p[j]);
     }
   }
   return arma::eig_sym(symMat);
 }

/***
 * Computes the asymptotic CDF function in the mixed case.
 */
// [[Rcpp::export]]
arma::vec HoeffIndMixedCdfRCPP(arma::vec x, arma::vec eigenP, double maxError) {
  AsymMixedCdfIntegrandEvaluator amcie(eigenP);
  arma::vec cdfVals(x.size());
  for (unsigned int i = 0; i < x.size(); i++) {
    cdfVals[i] = boundInZeroOne(
      numericalCfInversion(amcie, x[i], 10.0, maxError, 12));
  }
  return cdfVals;
}
