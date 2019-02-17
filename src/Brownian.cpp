#include "Brownian.hpp"

Brownian::Brownian(double h, double T)
    : Path1D(h, T), h(h), T(T), G(0, 1)
{

}

std::vector<double> & Brownian::operator()(std::mt19937_64 & gen) {
    for (unsigned i = 1; i < values.size(); ++i)
        values[i] = values[i - 1] + std::sqrt(h) * G(gen);
    return values;
}