/*
 * Aggregates the features of a stochastic homogeneous diffusion.
 * Inspired by Vincent Lemaire's implementation, adapted for own enjoyment.
*/

#include <functional>

template <typename TState = double, typename TSigma = double>
class SDE {

public:
    SDE(TState const & init_value) : init_value(init_value) { }
    
    virtual TState b(TState const &) = 0;
    virtual TSigma sigma(TState const &) = 0;

protected:
    TState init_value;
};
