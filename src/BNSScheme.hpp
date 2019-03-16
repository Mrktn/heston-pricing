/*
 * The scheme associated with the BNS model presented by Pagès & Panloup.
 * The scheme is for (X, v) in this order.
*/

#include <random>
#include "BNSSDE.hpp"
#include "TSSFactory.hpp"

class BNSScheme {

public:
    // Cast operator
    operator std::array<double, 2>() const { return current_state; }

    BNSScheme(BNSSDE const & sde, std::function<double(unsigned)> const & gamma) : sde(sde), tss(sde.alpha, sde.c, sde.lambda), current_state(sde.init_value), gamma(gamma), n(0), G(0, 1) {}

    std::array<double, 2> operator()(std::mt19937_64 & gen) {
        double g = gamma(n + 1);
        double Z[] = {G(gen), G(gen)};
        double jump = tss(gen, g);
        double vt = current_state[1];

        current_state = std::array<double, 2> {
            current_state[0] + (sde.r - 0.5 * vt)*g + sde.rho * std::sqrt(vt * g) * Z[1] + std::sqrt((1 - sde.rho * sde.rho) * vt * g) * Z[0] + sde.phi * jump,
            std::abs(vt + sde.k * g * (sde.theta - vt) + sde.dzeta * std::sqrt(current_state[1] * g) * Z[1] + jump)
        };

        n += 1;
        return current_state;
    }

protected:
    BNSSDE sde;
    TSSFactory tss;
    std::array<double, 2> current_state;
    std::function<double(unsigned)> gamma;
    unsigned n;
    std::normal_distribution<> G;
};
