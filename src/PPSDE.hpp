/*
 * The Pagès-Panloup SDE yields a stationary process suitable for derivatives pricing
 * under Heston assumptions and decreasing Euler time steps.
*/

#include "SDE.hpp"

class PPSDE : public SDE<std::array<double, 2>, std::array<double, 2>> {

    // The associated scheme should be able to access protected fields
    friend class PPScheme;

public:
    PPSDE(std::array<double, 2> const & initial, double k, double theta, double dzeta) : SDE(initial), k(k), theta(theta), dzeta(dzeta) { }

    std::array<double, 2> b(std::array<double, 2> const & X) {
        return std::array<double, 2> { -X[0], k * (theta - X[1])};
    }

    std::array<double, 2> sigma(std::array<double, 2> const & X) {
        return std::array<double, 2> {std::sqrt(X[1]), dzeta * std::sqrt(X[1])};
    }

protected:
    double k, theta, dzeta;
};
