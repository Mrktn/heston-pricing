#include <vector>

/*
 * Generic class which represents 1-dimensional discrete paths of stochastic processes.
*/
class Path1D {

public:
    Path1D(double h, double T) : times(std::ceil(T / h)), values(std::ceil(T / h), 0.0) {
        for(unsigned k = 1; k < times.size(); ++k)
            times[k] = k * h * T;
    };

    // Returns the values of a newly generated path
    virtual std::vector<double> & operator()(std::mt19937_64 & gen) = 0;

    std::vector<double> times;
    std::vector<double> values;
};
