Projet C++, Probabilités Numériques pour la Finance
Projet 2.3 - Modèles a volatilité stationnaire I
Antoine Balestrat et Enoch Yang, promo 2018-2019

===================
=== HOW TO RUN ====
===================

mkdir bin
mkdir obj
make
./bin/heston-pricing [--flags, -h for help]

============================
=== MAIN FUNCTIONALITIES ===
============================

Main functionalities (flags), all will be called using default numeric values:

1. --type asian --model heston
2. --type asian --model bns
3. --type european --model heston
4. --type european --model heston --imply

===========================================
=== DIRTREE & DESCRIPTION OF MAIN FILES ===
===========================================

benchmarks/                   -- benchmark plots
	   test_iters.py
	   test_maturity.py
img/                          -- some plots
src/
    AsianCall.hpp             -- asian call pricing using heston model
    BNSAsianCall.hpp          -- asian call pricing using BNS model
    BNSSDE.hpp                -- class for the BNS SDE
    BNSScheme.hpp             -- discretization scheme of the BNS SDE
    Brownian.hpp              -- d-dimensional Brownian motion
    Brownian1D.hpp            -- 1 dimensional Brownian motion
    EuropeanCall.hpp          -- european call pricing using heston model
    PPSDE.hpp                 -- class for the SDE defined in Pages & Panloup
    PPScheme.hpp              -- discretization scheme of the above SDE
    SDE.hpp                   -- template for SDE's
    TSSFactory.hpp            -- temperate state process generation (for the BNS model)
    Trajectory.hpp            -- container for d-dimensional trajectories
    Trajectory1D.hpp          -- container for 1-dimensional trajectories
    black_scholes.hpp         -- B-S european call
    bs_formulas.hpp           -- B-S formulas
    commandline.hpp           -- header for the command line interface (CLI)
    commandline.cpp           -- CLI
    main.cpp
    vol_implier.hpp           -- class to calculate the implied volatility

Makefile                      -- run make to compile the code
README.txt
hist.py
