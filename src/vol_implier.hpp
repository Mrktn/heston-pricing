/*
 * Functions to generate the implied volatility by an european call price
 */

#include "black_scholes.hpp"
#include <cmath>
#include <boost/math/tools/roots.hpp>

template <class T>
struct bscall_vol_functor
{
    bscall_vol_functor(BlackScholesCall<T>& bs, T price) :
	     bs(bs), price(price) { };
    
    std::tuple<T, T, T> operator()(T const& x)
    {
        T fx = bs.option_price(x) - price;
	T dx = bs.option_vega(x);
	T d2x = bs.option_vomma(x);
    	return std::make_tuple(fx, dx, d2x);
    }

private:
    BlackScholesCall<T> bs;
    T price;
};

template <class T>
class VolImplier {

public:
    VolImplier(T S, T K, T r, T M) :
	    bs(BlackScholesCall<T>(S, K, r, M))
    {
	// set initial sigma guess at inflexion point
	// Manaster-Koehler 1982
	initial_sigma = sqrt(2 * abs(log(S*exp(r*M)/K)) / M);	
    };

    T implied_vol(const T price)
    {
	
        using namespace boost::math::tools;
	T min = 0.0, max = 1.0;
	const int digits = std::numeric_limits<double>::digits;
	int get_digits = static_cast<int>(digits * 0.5);

	boost::uintmax_t maxit = 20;
	T result = halley_iterate(bscall_vol_functor<T>(bs, price), initial_sigma,
			min, max, get_digits, maxit);
	return result;
    }

private:
    BlackScholesCall<T> bs;
    T initial_sigma;
};
