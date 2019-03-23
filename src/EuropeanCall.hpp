/*
 * Approximates the price of a vanilla European Call option.
*/

#include <cassert>
#include <array>
#include <functional>
#include "PPScheme.hpp"

class EuropeanCall {

public:
    EuropeanCall(double s0, double r, double T, double K, double rho, double k, double theta, double dzeta) : s0(s0), r(r), T(T), K(K), rho(rho), k(k), theta(theta), dzeta(dzeta) {
        // Asserts that the Berkaoui et al. condition is satisfied
        //assert(2*k*theta / (dzeta*dzeta) > 1 + 2 * std::sqrt(6) / dzeta);
    }

    // TODO: successive calls of an iter() method ?
    // n_iter
    double simulate(unsigned n_iter, std::function<double(unsigned)> const & gamma, std::function<double(unsigned)> const & eta, std::mt19937_64 & gen) {
        double Gamma_running = 0.0, Gamma_ref = 0.0;
        unsigned curr = 0;

        std::array<double, 2> const init = {0.0, theta};
        PPSDE sde(init, k, theta, dzeta);
        PPScheme scheme(&sde, gamma);

        const double kappa      = r - rho * k * theta / dzeta;
        const double beta       = rho * k / dzeta - 0.5;
        const double mu         = rho / dzeta;
        const double upsilon    = std::sqrt(1 - rho * rho);

        double prev_Gamma = 0.0;
        std::array<double, 2>  prev_X = init;
        std::array<double, 2>  pprev_X;
        double prev_M = 0.0, prev_V = 0.0;
        double pprev_M = 0.0;
	//double pprev_V = 0.0;
        double prev_psi = s0;
        double curr_nu = 0.0;
        double H = 0.0;
        double psi = 0.0;

        std::vector<double> Xi = {1.0};
        
        double last_increment_psi = 0.0;
	//double last_increment_t = 0.0;

        for (unsigned k_iter = 0; k_iter < n_iter; ++k_iter) {
            do {
                curr++;
                Gamma_running += gamma(curr);

                std::array<double, 2> const X = scheme(gen);
                double deltat = Gamma_running - prev_Gamma;
                double M = X[0];
                double V = prev_V + deltat * prev_X[1];
                double increment_psi = std::exp((kappa * deltat) + (beta * prev_X[1] * deltat) + (mu * (X[1] - prev_X[1])) + (upsilon * (M - prev_M)));
                psi = prev_psi * increment_psi;
                //last_increment_t = deltat;
                last_increment_psi = increment_psi;

                if (curr <= n_iter)      
                {    
                    Xi.push_back(std::exp(kappa * Gamma_running + beta * V + upsilon * M));
                };

                prev_Gamma = Gamma_running;
                pprev_M = prev_M;
                //pprev_V = prev_V;
                pprev_X = prev_X;
                prev_M = M;
                prev_V = V;
                prev_X = X;
                prev_psi = psi;

            } while (Gamma_running - Gamma_ref <= T);

            double amend_deltat = Gamma_running - (T + Gamma_ref);
            double amended_truePsi = psi / last_increment_psi;
            amended_truePsi *= std::exp((kappa * amend_deltat) + (beta * pprev_X[1] * amend_deltat) + (mu * (prev_X[1] - pprev_X[1])) + (upsilon * (prev_M - pprev_M)));

            amended_truePsi /= Xi[k_iter];
            //double F = std::exp(-r*T) * std::max(amended_truePsi - K, 0.0);
            double F = s0 - std::exp(-r * T)*K + std::exp(-r * T) * std::max(K - amended_truePsi, 0.0); // use call put parity to reduce variance

            H += eta(k_iter + 1);
            curr_nu = curr_nu + (eta(k_iter + 1) / H) * (F - curr_nu);

            Gamma_ref += gamma(k_iter + 1);
        }

        return curr_nu;
    }


protected:
    // Usual market variables
    double s0, r, T;

    // Strike
    double K;

    // Correlation between the stochastic drivers of the underlying and its spot variance
    double rho;

    // Heston parameters
    double k, theta, dzeta;

    //std::normal_distribution<> G;
};
