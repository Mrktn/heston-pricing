/*
 * This is the Pagès-Panloup (BNS-flavored) SDE with jumps as given in the reference paper.
 *
 * This is the SDE satisfied by (X, v)
 * dXt = (r - vt/2)dt + <(rho sqrt(vt), sqrt(1 - rho^2) sqrt(vt)), (dW^1_t, dW^2_t)> + phi dZt
 * dvt = k(theta - vt)dt + dzeta sqrt(vt) dW^2_t + dZt
*/

#include <array>

#include "SDE.hpp"

class BNSSDE : public SDE<std::array<double, 2>, std::array<std::array<double, 2>, 2>> {

    // The associated scheme should be able to access protected fields
    friend class BNSScheme;

public:
    BNSSDE(std::array<double, 2> const & initial, double k, double theta, double dzeta, double r, double rho, double alpha, double c, double lambda, double phi) : SDE(initial), k(k), theta(theta), dzeta(dzeta), r(r), rho(rho), alpha(alpha), c(c), lambda(lambda), phi(phi) { }

    std::array<double, 2> b(std::array<double, 2> const & X) {
        return std::array<double, 2> {r - X[1] / 2, k * (theta - X[1])};
    }

    std::array<std::array<double, 2>, 2> sigma(std::array<double, 2> const & X) {
        return std::array<std::array<double, 2>, 2> {
            std::array<double, 2> {rho * std::sqrt(X[1]), std::sqrt((1 - rho * rho)) * X[1]},
            std::array<double, 2> {0, dzeta * std::sqrt(X[1])}
        };
    }

    std::array<double, 2> kappa(std::array<double, 2> const & X) {
        return std::array<double, 2> {phi, 1};
    }

protected:
    double k, theta, dzeta;
    double r, rho;
    double alpha, c, lambda;
    double phi;
};
