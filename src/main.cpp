#include <iostream>
#include <random>
#include <algorithm>

#include "BNSAsianCall.hpp"

int main(int argc, char **argv) 
{
    std::random_device rd;
    std::mt19937_64 gen(rd());

    //AsianCall asset(50, 0.05, 1., 44, 0.5, 2, 0.01, 0.1);
    // BNSAsianCall(s0, r, T, K, rho, k, theta, dzeta, alpha, c, lambda, phi)
    BNSAsianCall asset(50, 0.05, 1, 44, 0.0, 1, 0.00, 0, 0.5, 0.01, 1, -1.0); 
    std::function<double(unsigned)> const gamma = [](unsigned n) {return 0.01*std::pow(n, - 1.0 / 3.0);};
    std::function<double(unsigned)> const eta   = [](unsigned n) {return std::pow(n, - 1.0 / 3.0);};

    std::cout << asset.simulate(5*100000, gamma, eta, gen) << std::endl;

    return 0;
}
