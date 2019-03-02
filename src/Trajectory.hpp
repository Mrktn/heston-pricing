#include <array>
#include <vector>
#include <functional>

/*
 * Generic class which represents d-dimensional discrete paths of stochastic processes, with potentially decreasing time steps.
*/
template <unsigned Dimension>
class Trajectory {

public:
    // Constructor for non uniform time steps
    Trajectory(const std::function<double(unsigned)> & gamma, unsigned n) : d(Dimension), n(n), gamma(gamma), times(n), values(n, std::array<double, Dimension>()) {
        double S = 0;
        for(unsigned k = 1; k < n; ++k) {
            S = S + gamma(k);
            times[k] = S;
        }
    };

    // Returns the values of a newly generated path
    virtual std::vector<std::array<double, Dimension>> & operator()(std::mt19937_64 & gen) = 0;


    // Dimension
    unsigned d;
    // How many points in the trajectory ?
    unsigned n;

    // In case of non uniform time steps, this is the (lowercase) gamma sequence : a function which maps n (unsigned) to gamma_n (double).
    std::function<double(unsigned)> gamma;

    // A trajectory is basically a vector of times and a vector of associated values
    std::vector<double> times;
    std::vector<std::array<double, Dimension>> values;
};
