/*
 * Generates paths of a 1D Brownian motion with unitary volatility, over [0 ; T].
*/

#include <vector>
#include <random>

#include "Trajectory1D.hpp"

class Brownian1D : public Trajectory1D {

public:
    // Constructor for uniform time steps
    Brownian1D(double h, unsigned n) : Trajectory1D(h, n), G(0, 1) { };

    // Constructor for non uniform time steps
    Brownian1D(std::function<double(unsigned)> const & gamma, unsigned n) : Trajectory1D(gamma, n), G(0, 1) { };

    // Returns the values of a newly generated path
    std::vector<double> & operator()(std::mt19937_64 & gen) {
        for (unsigned i = 1; i < values.size(); ++i)
            values[i] = values[i - 1] + std::sqrt(times[i] - times[i - 1]) * G(gen);
        return values;
    };

protected:
    std::normal_distribution<> G;
};
