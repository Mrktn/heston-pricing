#include <iostream>
#include <random>
#include <algorithm>

#include "Brownian.hpp"
#include "PPScheme.hpp"

int main(int argc, char **argv) 
{
    std::random_device rd;
    std::mt19937_64 gen(rd());

    PPSDE sde(std::array<double, 2> {1.0, 0.0}, 2.0, 0.01, 0.1);

    PPScheme s(sde, [](unsigned n) { return std::pow(n, - 1.0 / 3.0); }, 100);

    return 0;
}
