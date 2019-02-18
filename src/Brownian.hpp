/*
 * Generates paths of a 1D Brownian motion with unitary volatility, over [0 ; T].
*/

#include <vector>
#include <random>

#include "Trajectory.hpp"

template <unsigned Dimension>
class Brownian : public Trajectory<Dimension> {

public:
    // Constructor for uniform time steps
    Brownian(double h, unsigned n) : Trajectory<Dimension>(h, n), G(0, 1) { };

    // Constructor for non uniform time steps
    Brownian(const std::function<double(unsigned)> & gamma, unsigned n) : Trajectory<Dimension>(gamma, n), G(0, 1) { };

    // Returns the values of a newly generated path
    std::vector<std::array<double, Dimension>> & operator()(std::mt19937_64 & gen) {
        for (unsigned i = 1; i < this->values.size(); ++i) {
            for (unsigned k = 1; k < Dimension; ++k) {
                this->values[i][k] = this->values[i - 1][k] + std::sqrt(this->times[i] - this->times[i - 1]) * G(gen);
            }
        }
        return this->values;
    };

protected:
    std::normal_distribution<> G;
};
