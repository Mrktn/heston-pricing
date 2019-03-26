/*
 * The Berkaoui et al scheme as presented by Pagès & Panloup.
 * The scheme is for (M, v) in this order.
*/

#ifndef PPSCHEME_HPP_INCLUDED
#define PPSCHEME_HPP_INCLUDED

#include <random>
#include "PPSDE.hpp"

class PPScheme {

public:
    // Cast operator
    operator std::array<double, 2>() const { return current_state; }

    PPScheme(PPSDE * sde, std::function<double(unsigned)> const & gamma) : sde(sde), current_state(sde->init_value), gamma(gamma), n(0), G(0, 1) {}

    std::array<double, 2> operator()(std::mt19937_64 & gen) {
        double g = gamma(n + 1);
        double Z[] = {G(gen), G(gen)};

        current_state = std::array<double, 2> {
            current_state[0] + std::sqrt(current_state[1] * g) * Z[0],
            std::abs(current_state[1] + sde->k * g * (sde->theta - current_state[1]) + sde->dzeta * std::sqrt(current_state[1] * g) * Z[1])
        };

        n += 1;
        return current_state;
    }

protected:
    PPSDE *sde;
    std::array<double, 2> current_state;
    std::function<double(unsigned)> gamma;
    unsigned n;
    std::normal_distribution<> G;
};

#endif