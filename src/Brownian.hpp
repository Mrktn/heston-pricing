/*
 * Generates paths of a Brownian motion with unitary volatility, over [0 ; T].
*/

#include <vector>
#include <random>

#include "Path1D.hpp"

class Brownian : public Path1D {

public:
    Brownian(double h = 0.5, double T = 1);

    // Returns the values of a newly generated path
    std::vector<double> & operator()(std::mt19937_64 & gen);

protected:
    double h, T;
    std::normal_distribution<> G;
};
