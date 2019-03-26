/*
 * Aggregates the features of a homogeneous Lévy-driven SDE.
 * Inspired by Vincent Lemaire's implementation, adapted for own enjoyment.
 *
 * dXt = b(Xt)dt + sigma(Xt)dWt + k(Xt)dZt
 * where Z is an integrable, purely discontinuous, R-valued Lévy process
*/

#ifndef SDE_HPP_INCLUDED
#define SDE_HPP_INCLUDED

#include <functional>

template <typename TState = double, typename TSigma = double>
class SDE {

public:
    SDE(TState const & init_value) : init_value(init_value) { }
    
    virtual TState b(TState const &) = 0;
    virtual TSigma sigma(TState const &) = 0;
    virtual TState kappa(TState const &) = 0;

protected:
    TState init_value;
};

#endif