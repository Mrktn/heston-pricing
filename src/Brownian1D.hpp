/*
 * Generates paths of a 1D Brownian motion with unitary volatility, over [0 ; T].
*/

#include <vector>
#include <random>

#include "Trajectory1D.hpp"

class Brownian1D : public Trajectory1D {

public:
    // Constructor for uniform time steps
    Brownian1D(double h, unsigned n);

    // Constructor for non uniform time steps
    Brownian1D(const std::function<double(unsigned)> & gamma, unsigned n);

    // Returns the values of a newly generated path
    std::vector<double> & operator()(std::mt19937_64 & gen);

protected:
    std::normal_distribution<> G;
};
