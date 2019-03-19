/*
 * Generates independent increments of a spectrally positive tempered stable process
 * R. Kawai, H. Masuda, On Simulation of Tempered Stable Random Variates
*/

#include <cassert>
#include <cmath>
#include <random>

class TSSFactory {

public:

    TSSFactory(double alpha, double a, double b) : alpha(alpha), a(a), b(b), unif_pi(- M_PI / 2.0, M_PI / 2), unif_n(0.0, 1.0), expon(1.0) {
        assert((0 < alpha && 1 > alpha) && "Can't guarantee the rejection method sampling for alpha not in (0, 1) !");
    }

    // Rejection method for the computation of the increment
    double operator()(std::mt19937_64 & gen, double deltat) {
        bool done = false;
        double ret;

        while (!done) {
            double U = unif_n(gen);
            double V = S(deltat, gen);

            if (U <= std::exp(- b * V)) {
                done = true;
                ret = V;
            }
        }

        return ret;
    }

protected:
    double alpha, a, b;
    std::uniform_real_distribution<double> unif_pi, unif_n;
    std::exponential_distribution<double> expon;

private:
    // S(alpha, a * deltat)
    double S(double deltat, std::mt19937_64 & gen) {
        double U = unif_pi(gen), E = expon(gen);
        double aprime = a*deltat;
        double Theta = std::atan(std::tan(M_PI * alpha / 2.0));

        double X = std::pow(- aprime * std::tgamma(- alpha) * std::cos(M_PI * alpha / 2.0), 1.0 / alpha);
        double Y = std::sin(alpha * U + Theta) / std::pow(std::cos(U) * std::cos(Theta), 1.0 / alpha);
        double Z = std::pow(std::cos((1.0 - alpha) * U - Theta) / E, (1.0 - alpha) / alpha);

        return X * Y * Z;
    }
};
