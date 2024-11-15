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

#ifndef PseudotimeDE_IntegrandEvaluator
#define PseudotimeDE_IntegrandEvaluator

#include <complex>

/***
 * A class whose only purpose is to be subclassed. Provides an interface
 * to get an integrand (taking two parameters) up to some maximum error.
 */
class IntegrandEvaluator { 
public:
  virtual std::complex<double> integrand(double x, double t, double maxError) = 0;
};

#endif
