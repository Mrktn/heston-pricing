/*
 * The Berkaoui et al scheme as presented by Pagès & Panloup.
 * The scheme is for (y, v) in this order.
*/

// TODO: make it a subclass of a more general "Scheme" class

#include "PPSDE.hpp"

class PPScheme {

public:
    // Cast operator
    operator std::array<double, 2>() const { return current_state; }

    PPScheme(PPSDE const & sde, std::function<double(unsigned)> const & gamma) : sde(sde), current_state(sde.init_value), gamma(gamma), n(0) {}

    ~PPScheme() {std::cout << "oskour depuis scheme" << std::endl;}

    std::array<double, 2> operator()(std::array<double, 2> const & Z) {
        double g = gamma(n + 1);

        current_state = std::array<double, 2> {
            current_state[0] - g*current_state[0] + std::sqrt(current_state[1] * g) * Z[0],
            std::abs(current_state[1] + sde.k * g * (sde.theta - current_state[1]) + sde.dzeta * std::sqrt(current_state[1] * g) * Z[1])
        };

        n += 1;
        return current_state;
    }

protected:
    PPSDE sde;
    std::array<double, 2> current_state;
    std::function<double(unsigned)> gamma;
    unsigned n;
};
