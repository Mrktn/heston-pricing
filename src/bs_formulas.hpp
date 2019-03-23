/*
 * Small library of functions to calculate the Black Scholes formula and some Greeks
 * TODO maybe make this a template
 */ 

#ifndef __BS_FORMULAS_H
#define __BS_FORMULAS_H

#define __USE_MATH_DEFINES

#include <iostream>
#include <cmath>

double normal_pdf(const double x)
{
    return (1.0 / (pow(2*M_PI, 0.5)))*exp(-0.5*x*x); 
}

double normal_cdf(const double x)
{
    return 0.5 * erfc(-x * M_SQRT1_2);
}

double d_j(const int j, const double S, const double K, const double r,
		const double sigma, const double T)
{
    return (log(S/K) + (r + (pow(-1,j-1))*0.5*sigma*sigma)*T)
	    /(sigma*(sqrt(T)));
}

double call_price(const double S, const double K, const double r,
		const double sigma, const double T)
{
  return S * normal_cdf(d_j(1, S, K, r, sigma, T))
	  -K*exp(-r*T) * normal_cdf(d_j(2, S, K, r, sigma, T));
}

double call_vega(const double S, const double K, const double r,
		const double sigma, const double T)
{
  return S * normal_pdf(d_j(1, S, K, r, sigma, T))*sqrt(T);
  // return K * exp(-r*T) * normal_pdf(d_j(2, S, K, r, sigma, T)*pow(T, 0.5);
}

double call_vomma(const double S, const double K, const double r,
		const double sigma, const double T)
{
    return call_vega(S, K, r, sigma, T) * d_j(1, S, K, r, sigma, T) *
	    d_j(2, S, K, r, sigma, T) / sigma;
}

#endif
