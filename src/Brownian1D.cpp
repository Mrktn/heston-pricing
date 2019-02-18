#include "Brownian1D.hpp"

Brownian1D::Brownian1D(double h, unsigned n)
    : Trajectory1D(h, n), G(0, 1)
{

}

Brownian1D::Brownian1D(const std::function<double(unsigned)> & gamma, unsigned n)
    : Trajectory1D(gamma, n), G(0, 1)
{

}

std::vector<double> & Brownian1D::operator()(std::mt19937_64 & gen) {
    for (unsigned i = 1; i < values.size(); ++i)
        values[i] = values[i - 1] + std::sqrt(times[i] - times[i-1]) * G(gen);
    return values;
}