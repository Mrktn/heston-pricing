#ifndef __BLACK_SCHOLES_H
#define __BLACK_SCHOLES_H

#include "bs_formulas.hpp"

template <class T>
struct BlackScholesCall 
{

    BlackScholesCall(T S, T K, T r, T M) :
	    S(S), K(K), r(r), M(M) { };

    T option_price(T const& sigma) const {
    	return call_price(S, K, r, sigma, M);
    }

    T option_vega(T const& sigma) const {
    	return call_vega(S, K, r, sigma, M);
    }

    T option_vomma(T const& sigma) const {
    	return call_vomma(S, K, r, sigma, M);
    }

private:
    T S, K, r, M;
};

#endif
