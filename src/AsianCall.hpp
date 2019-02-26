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
        //current_n = 0;
    }



    // Computes F over the X = (y, v) trajectory
    // Basically : (exp(-rT) / T * integral - K)^+
    /*double F(std::vector<double> const & times, std::vector<std::array<double, 2>> const & X) {
        //return std::max(0.0, std::exp(-r*T) * integral(times, X) / T - K);
        return 0.0;
    }*/

    // TODO: successive calls of an iter() method ?
    // n_iter
    double simulate(unsigned n_iter, std::function<double(unsigned)> const & gamma, std::function<double(unsigned)> const & eta, std::mt19937_64 & gen) {
        unsigned curr_nu_n = 0;
        double curr_nu = 0.0;
        double H = 0;

        N_total_iters = n_iter;

        // CHECKME: does initializing v0 to theta make a difference ? Intuitively, should speed up the process even by a little...
        std::array<double, 2> const init = {0.5, theta};

        saved_elements = new double[n_iter+1];
        saved_elements[0] = s0;
        last_undl = s0;

        PPSDE sde(init, k, theta, dzeta);
        PPScheme scheme(sde, gamma);

        unsigned k_iter = 0;
        double Gamma_act = 0, Gamma_prev = 0;

        std::vector<std::array<double, 2>> X;
        std::vector<double> times;

        X.push_back(init);
        times.push_back(0.0);

        unsigned start_integral = 0;
        double prev_integral = 0.0;
        double retrieve = 0.0;



        ////////////////////////////////////
            std::vector<double> loc_times;
            std::vector<std::array<double, 2>> loc_X;
            loc_X.push_back(init);
            loc_times.push_back(0.0);
            double loc_S = 0.0;
            for(unsigned k = 1; k < n_iter; ++k) {
                loc_S += gamma(k);
                loc_times.push_back(loc_S);
                loc_X.push_back(scheme(std::array<double, 2> {G(gen), G(gen)}));
                
            }
            integral(loc_times, loc_X, 0);
            return 0.0;

            ////////////////////////////////////

        while (curr_nu_n < n_iter) {
            //std::cout << std::endl;
            //std::cout << "Pour curr_nu = " << curr_nu_n << std::endl;

            //std::cout << "Gamma actuel = " << Gamma_act << std::endl;
            //std::cout << "Gamma previous = " << Gamma_prev << std::endl;
            // First : iterate X's scheme far enough
            while (Gamma_act - Gamma_prev <= T) {
                X.push_back(scheme(std::array<double, 2> {G(gen), G(gen)}));
                Gamma_act += gamma(k_iter + 1);
                times.push_back(Gamma_act);
                ++k_iter;
            }

            Gamma_prev += gamma(curr_nu_n + 1);

            // J'ai [X_G_0, X_G_1, ..., X_G_N(0, T)] dans X
            auto X_bckp = X.back(); auto t_bckp = times.back();
            X.pop_back();
            times.pop_back();


            /*std::cout << "Calling integral with" << std::endl;
            for(size_t k = 0; k < times.size(); ++k)
                std::cout << times[k] << std::endl;*/

            double integr = integral(times, X, start_integral);
            double actual = prev_integral + integr - retrieve;
            /*std::cout << "Actual : " << actual << std::endl;
            std::cout << "integr : " << integr << std::endl;
            std::cout << "Prev : " << prev_integral << std::endl;
            std::cout << "Retrieve : " << retrieve << std::endl;*/
            double F_tmp = std::max(0.0, std::exp(-r * T) * actual / T - K);
            //std::cout << "F_tmp : " << F_tmp << std::endl;
            prev_integral = actual;
            retrieve = saved_elements[curr_nu_n];
            start_integral += times.size()-1;

            X.clear(); X.push_back(X_bckp);
            times.clear(); times.push_back(t_bckp);

            curr_nu_n += 1;
            double et = eta(curr_nu_n);
            H += et;

            //std::cout << "Allez zou, on rajoute " << et << " et on devient H = " << H << std::endl;
            curr_nu = curr_nu + (et / H) * (F_tmp - curr_nu);
            //std::cout << "curr_nu est maintenant " << curr_nu << std::endl;
            // k_iter vaut N(0, T)
        }

        /*std::cout << "Mon underlying :" << std::endl;
        for(double *s = underlying; s < underlying + n_iter; s += 1)
            std::cout << *s << std::endl;*/

        delete[] saved_elements;
        return curr_nu;
    }


private:

    unsigned N_total_iters;

    double *saved_elements;
    double last_undl;

    // Integral computed from start
    double integral(std::vector<double> const & times, std::vector<std::array<double, 2>> const & X, size_t start) {
        std::cout << "Called integral" << std::endl;
        // M[i] is an approximation of the special (stationary) transform of y = X[0] at time i as discussed in Pagès & Panloup
        std::vector<double> M(times.size(), 0.0);

        // V[i] is an approximation of the integral of v = X[1] till te ith time
        std::vector<double> V(times.size(), 0.0);

        // psi[i] is the spot price at the ith time
        //std::cout << "Attention, j'utilise underlying au temps " << start << " où il vaut " << last_undl << std::endl;
        std::vector<double> psi(times.size(), last_undl);

        // ret is an discrete approximation of the payoff of the asian call
        double ret = 0;

        //std::cout << "start = " << start << std::endl;
        // TODO: merge these loops once everything works fine

        for (unsigned j = 1; j < times.size(); ++j) {
            M[j] = M[j - 1] + (X[j][0] - X[j - 1][0]) + X[j - 1][0] * (times[j] - times[j - 1]);
            V[j] = V[j - 1] + (times[j] - times[j - 1]) * X[j - 1][1];
        }

        for (unsigned j = 1; j < times.size(); ++j) {
            double deltat = times[j] - times[j - 1];
            double deltaV = V[j] - V[j - 1];
            
            psi[j] = psi[j - 1] * std::exp(r * deltat - deltaV + (rho / dzeta) * (X[j][1] - X[j - 1][1] - k * theta * deltat + k * deltaV) + std::sqrt(1 - rho * rho) * (M[j] - M[j - 1]));
            std::cout << times[j] << " " << psi[j] << std::endl;

            //std::cout << "underlying est calculé au temps " << start + j << " et vaut " << psi[j] << std::endl;

            if(j == times.size() - 1) {
                //std::cout << "Le dernier temps calculé est donc t = " << start + j << " et s = " << psi[j] << std::endl;
                last_undl = psi[j];
            }
            //psi[j] = psi[j - 1] * std::exp(r * deltat - deltaV + (rho / dzeta) * (X[j][1] - X[j - 1][1] - k * theta * deltat + k * deltaV) + std::sqrt(1 - rho * rho) * (M[j] - M[j - 1]));

            

            //std::cout << times[j] << " " << psi[j] << std::endl;
        }

        // The superfluous term is meant to be substracted in order to compute recursively the integral over the very next path
        //superfluous = psi[0] * (times[1] - times[0]);

        for (unsigned k = 0; k < times.size()-1; ++k) {
            double add = psi[k] * (times[k + 1] - times[k]);
            ret += add;
            // The integral function is responsible for populating the underlying trajectory array (FIXME: can it be done elsewhere ?)
            if (start + k < N_total_iters) {
                //std::cout << "Setting " << start + k << std::endl;
                saved_elements[start + k] = add;
            }
        }

        if (start + times.size()-1 < N_total_iters) {
                //std::cout << "Setting " << start + times.size()-1 << std::endl;
                saved_elements[start + times.size()-1] = psi[times.size()-1] * (times[times.size()-1 + 1] - times[times.size()-1]);
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
