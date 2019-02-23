/*
 * Approximates the price of an asian call.
 * TODO: the prices for several strikes can be computed at the same time in a single simulation
 * TODO: should 1) inherit common features from a superclass (Derivative ? PriceablePayoff ?) 
 *              2) befriend a HestonUnderlying class  
*/

#include <cassert>

class AsianCall {

public:
    AsianCall(double s0, double r, double T, double K, double rho, double k, double theta, double dzeta) : s0(s0), r(r), T(T), K(K), rho(rho), k(k), theta(theta), dzeta(dzeta), G(0, 1) {
        // Asserts that the Berkaoui et al. condition is satisfied
        //assert(2*k*theta / (dzeta*dzeta) > 1 + 2 * std::sqrt(6) / dzeta);
        current_n = 0;
    }



    // Computes F over the X = (y, v) trajectory
    // Basically : (exp(-rT) / T * integral - K)^+
    double F(std::vector<double> const & times, std::vector<std::array<double, 2>> const & X) {
        return std::max(0.0, std::exp(-r*T) * integral(times, X) / T - K);
    }

    // TODO: successive calls of an iter() method ?
    double simulate(unsigned n, std::function<double(unsigned)> const & gamma, std::function<double(unsigned)> const & eta, std::mt19937_64 & gen) {
        // CHECKME: does initializing v0 to theta make a difference ? Intuitively, should speed up the process even by a little...
        std::array<double, 2> const init = {0.5, theta};

        PPSDE sde(init, k, theta, dzeta);
        PPScheme scheme(sde, gamma);

        unsigned k_iter = 0;
        double Gamma_act = 0, Gamma_prev = 0;

        std::vector<std::array<double, 2>> X;
        std::vector<double> times;
        
        X.push_back(init);
        times.push_back(0.0);

        while(Gamma_act - Gamma_prev <= T) {
            X.push_back(scheme(std::array<double, 2> {G(gen), G(gen)}));
            Gamma_act += gamma(k_iter + 1);
            times.push_back(Gamma_act - Gamma_prev);
            ++k_iter;
        }

        // J'ai [X_G_0, X_G_1, ..., X_G_N(0, T)] dans X
        X.pop_back();
        times.pop_back();
        unsigned N = k_iter - 1;

        return F(times, X);

        // k_iter vaut N(0, T)
    }


private:
    // What is the last nu approximation that we computed ?
    unsigned current_n;

    double integral(std::vector<double> const & times, std::vector<std::array<double, 2>> const & X) {
        // M[i] is an approximation of the special (stationary) transform of y = X[0] at time i as discussed in Pagès & Panloup
        std::vector<double> M(times.size(), 0.0);

        // V[i] is an approximation of the integral of v = X[1] till te ith time
        std::vector<double> V(times.size(), 0.0);

        // psi[i] is the spot price at the ith time
        std::vector<double> psi(times.size(), s0);

        // ret is an discrete approximation of the payoff of the asian call
        double ret = 0;

        // TODO: merge these loops once everything works fine

        for(unsigned j = 1; j < times.size(); ++j) {
            M[j] = M[j-1] + (X[j][0] - X[j-1][0]) + X[j-1][0] * (times[j] - times[j-1]);
            V[j] = V[j-1] + (times[j] - times[j-1]) * X[j-1][1];
        }

        for(unsigned j = 1; j < times.size(); ++j) {
            double deltat = times[j] - times[j-1];
            double deltaV = V[j] - V[j-1];
            psi[j] = psi[j-1] * std::exp(r * deltat - deltaV + (rho / dzeta) * (X[j][1] - X[j-1][1] - k * theta * deltat + k * deltaV) + std::sqrt(1 - rho*rho) * (M[j] - M[j-1]));
            std::cout << times[j] << " " << psi[j] << std::endl;
        }

        for(unsigned k = 0; k < times.size() - 1; ++k) {
            ret += psi[k] * (times[k+1] - times[k]);
        }

        return ret;
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

    std::normal_distribution<> G;
/*
    // Number of iterations of the simulation
    unsigned n;

    // The step function
    std::function<double(unsigned)> gamma;
    // The weighting function
    std::function<double(unsigned)> eta;
*/
};
