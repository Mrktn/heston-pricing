#include <iostream>
#include <random>
#include <algorithm>

#include "Brownian.hpp"
#include "PPScheme.hpp"
#include "common.hpp"
#include "NAsianCall.hpp"
#include "AsianCall.hpp"

int main(int argc, char **argv) 
{
    std::random_device rd;
    std::mt19937_64 gen(rd());

    NAsianCall asset(50, 0.05, 1., 44, 0.5, 2, 0.01, 0.1);
    std::function<double(unsigned)> const gamma = [](unsigned n) {return std::pow(n, - 1.0 / 3.0);};
    std::function<double(unsigned)> const eta   = [](unsigned n) {return std::pow(n, - 1.0 / 3.0);};

    std::cout << asset.simulate(5*10000, gamma, eta, gen) << std::endl;
    return 0;
}
