#include "Brownian.hpp"

Brownian::Brownian(double h, unsigned n)
    : Path1D(h, n), G(0, 1)
{

}

Brownian::Brownian(const std::function<double(unsigned)> & gamma, unsigned n)
    : Path1D(gamma, n), G(0, 1)
{

}

std::vector<double> & Brownian::operator()(std::mt19937_64 & gen) {
    for (unsigned i = 1; i < values.size(); ++i)
        values[i] = values[i - 1] + std::sqrt(times[i] - times[i-1]) * G(gen);
    return values;
}