#include <iostream>
#include <random>
#include <algorithm>

#include "Brownian.hpp"

int main(int argc, char **argv) 
{
    std::random_device rd;
    std::mt19937_64 gen(rd());

    Brownian<2> B(0.5, 20);

    return 0;
}
