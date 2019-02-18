/*
 * Generates paths of a Brownian motion with unitary volatility, over [0 ; T], with uniform time steps.
*/

#include <vector>
#include <random>

#include "Path1D.hpp"

class Brownian : public Path1D {

public:
    // Constructor for uniform time steps
    Brownian(double h, unsigned n);

    // Constructor for non uniform time steps
    Brownian(const std::function<double(unsigned)> & gamma, unsigned n);

    // Returns the values of a newly generated path
    std::vector<double> & operator()(std::mt19937_64 & gen);

protected:
    std::normal_distribution<> G;
};
