#include "commandline.hpp"

#include "EuropeanCall.hpp"
#include "AsianCall.hpp"
#include "BNSAsianCall.hpp"

using namespace boost::program_options;

void on_age(int age)
{
    std::cout << "On age: " << age << '\n';
}

bool params_fit_for_simu(variables_map & vm) {
    double gammaexp = vm["gammaexponent"].as<double>(), etaexp = vm["etaexponent"].as<double>();

    if (gammaexp >= 0 && etaexp >= 0) {
        std::cout << "Invalid exponents for gamma_n and/or eta_n (" << gammaexp << ", " << etaexp << ") !" << std::endl;
        return false;
    }

    return true;
}
int handle_commandline(int argc, char **argv) {
    try {
        options_description desc{"Options"};
        desc.add_options()
        ("help,h", "Shows this help")
        ("iters,i", value<int>()->default_value(5 * 1e4), "Number of iterations to perform")
        ("gammaexponent", value<double>()->default_value(-1.0 / 3.0, "-0.33"), "Exponent of the gamma series")
        ("etaexponent", value<double>()->default_value(-1.0 / 3.0, "-0.33"), "Exponent of the eta series")
        ("type,t", value<std::string>()->default_value("asian"), "Any of: european, asian")
        ("model,m", value<std::string>()->default_value("heston"), "Any of: heston, bns")
        ("S0", value<double>()->default_value(50, "50"), "Initial value of the asset")
        ("strike,K", value<double>()->default_value(44, "44"), "Strike of the contract")
        ("maturity,T", value<double>()->default_value(1, "1"), "Maturity of the contract (in years)")
        ("interest,r", value<double>()->default_value(0.05, "0.05"), "Interest rate")
        ("correlation", value<double>()->default_value(0.5, "0.5"), "Correlation spot / variance (rho)")
        ("elasticity", value<double>()->default_value(2, "2"), "Speed of the mean-reverting variance feature (k)")
        ("meanrevert", value<double>()->default_value(0.01, "0.01"), "Long-term variance mean reversion (theta)")
        ("volofvar", value<double>()->default_value(0.1, "0.01"), "Volatility of the instantaneous variance (dzeta)")
        ("jumpstability", value<double>()->default_value(0.5, "0.5"), "Stability order of the jump process (alpha)")
        ("jumpscale", value<double>()->default_value(0.01, "0.01"), "Scale of the jump process (c)")
        ("jumptemperedexp", value<double>()->default_value(1, "1"), "Exponent of the tempered stable process (lambda)")
        ("jumpspotrelative", value<double>()->default_value(-1, "-1"), "Intensity of the spot jump relative to var jump (phi)");

        // BNSAsianCall asset(50, 0.05, 1, 44, 0.0, 1, 0.00, 0, 0.5, 0.01, 1, -1.0);
        // AsianCall asset(50, 0.05, 1., 44, 0.5, 2, 0.01, 0.1);
        // AsianCall(double s0, double r, double T, double K, double rho, double k, double theta, double dzeta)

        // BNSAsianCall(double rho, double k, double theta, double dzeta, double alpha, double c, double lambda, double phi)
        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);
        notify(vm);

        if (vm.count("help")) {
            std::cout << desc << std::endl;
            return 0;
        }

        if (!params_fit_for_simu(vm))
            return 1;


        double const gammaexp = vm["gammaexponent"].as<double>(), etaexp = vm["etaexponent"].as<double>();
        std::function<double(unsigned)> const gamma = [gammaexp](unsigned n) {return 0.01 * std::pow(n, gammaexp);};
        std::function<double(unsigned)> const eta   = [etaexp](unsigned n) {return std::pow(n, etaexp);};

        std::random_device rd;
        std::mt19937_64 gen(rd());

        if (vm["type"].as<std::string>() == "asian" && vm["model"].as<std::string>() == "heston") {
            // Market parameters
            double S0 = vm["S0"].as<double>(), K = vm["strike"].as<double>(), T = vm["maturity"].as<double>(), r = vm["interest"].as<double>();
            // Model parameters
            double rho = vm["correlation"].as<double>(), k = vm["elasticity"].as<double>(), theta = vm["meanrevert"].as<double>(), dzeta = vm["volofvar"].as<double>();
            AsianCall asset(S0, r, T, K, rho, k, theta, dzeta);
            std::cout << asset.simulate(vm["iters"].as<int>(), gamma, eta, gen) << std::endl;
            return 0;
        }

        else if (vm["type"].as<std::string>() == "asian" && vm["model"].as<std::string>() == "bns") {
            // Market parameters
            double S0 = vm["S0"].as<double>(), K = vm["strike"].as<double>(), T = vm["maturity"].as<double>(), r = vm["interest"].as<double>();
            // Heston basic parameters
            double rho = vm["correlation"].as<double>(), k = vm["elasticity"].as<double>(), theta = vm["meanrevert"].as<double>(), dzeta = vm["volofvar"].as<double>();
            // TSS parameters
            double alpha = vm["jumpstability"].as<double>(), c = vm["jumpscale"].as<double>(), lambda = vm["jumptemperedexp"].as<double>(), phi = vm["jumpspotrelative"].as<double>();
            BNSAsianCall asset(S0, r, T, K, rho, k, theta, dzeta, alpha, c, lambda, phi);
            std::cout << asset.simulate(vm["iters"].as<int>(), gamma, eta, gen) << std::endl;
            return 0;
        }

        else if (vm["type"].as<std::string>() == "european" && vm["model"].as<std::string>() == "heston") {
            // Market parameters
            double S0 = vm["S0"].as<double>(), K = vm["strike"].as<double>(), T = vm["maturity"].as<double>(), r = vm["interest"].as<double>();
            // Model parameters
            double rho = vm["correlation"].as<double>(), k = vm["elasticity"].as<double>(), theta = vm["meanrevert"].as<double>(), dzeta = vm["volofvar"].as<double>();
            EuropeanCall asset(S0, r, T, K, rho, k, theta, dzeta);
            std::cout << asset.simulate(vm["iters"].as<int>(), gamma, eta, gen) << std::endl;
            return 0;
        }
        else {
            std::cout << "Invalid / unimplemented combination of payoff / model (" << vm["type"].as<std::string>() << " / " << vm["model"].as<std::string>() << ") !" << std::endl;
            return 1;
        }

    } catch (const error &ex) {
        std::cerr << ex.what() << '\n';
        return 256;
    }
}
