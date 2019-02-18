#include <iostream>
#include <random>
#include <algorithm>

#include "Brownian1D.hpp"

int main(int argc, char **argv) 
{
    std::random_device rd;
    std::mt19937_64 gen(rd());

    Brownian1D B(0.5, 10);
    B(gen);

    for(unsigned i = 0; i < B.times.size(); ++i)
        std::cout << B.times[i] << " " << B.values[i] << "\n";

    std::cout << "\n";

    const std::function<double(unsigned)> gamma = [](unsigned n) {return std::pow(1.0 / n, 1.0 / 3.0);};

    Brownian1D B2(gamma, 100);
    B2(gen);

    for(unsigned i = 0; i < B2.times.size(); ++i)
        std::cout << B2.times[i] << " " << B2.values[i] << "\n";

    return 0;
}
