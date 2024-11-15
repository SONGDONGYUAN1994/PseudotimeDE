/***
 * Copyright (C) 2016 Luca Weihs
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

#ifndef PseudotimeDE_AsymMixedCdfIntegrandEvaluator
#define PseudotimeDE_AsymMixedCdfIntegrandEvaluator

// #include "AsymCdfIntegrandEvaluator.h"
#include "RcppArmadillo.h"
#include "IntegrandEvaluator.h"

double hurwitzZeta(double exponent, double offset, double maxError);
std::complex<double> gridSum(std::complex<double> v, int sideLen);
std::complex<double> sinhProd(std::complex<double> v, int i);
std::complex<double> asymContCharFunction(double t, double maxError);
double aCoef(int k, int h, double maxError);
std::complex<double> tailSum(std::complex<double> v, int h, double maxError);

int piRemSign(double x);
int getSinhSign(double rate);

class AsymMixedCdfIntegrandEvaluator : public IntegrandEvaluator {
protected:
  arma::vec eigenP;

public:
  AsymMixedCdfIntegrandEvaluator(arma::vec eigP);
  std::complex<double> integrand(double x, double t, double maxError);
};

#endif
