#!/usr/bin/python
from subprocess import Popen, PIPE
import numpy as np
import matplotlib.pyplot as plt 
import time
plt.style.use('seaborn-whitegrid')

iters = range(1, 500)

x = []
yheston = []
ybns = []

for i in iters:
    # Asian Heston
    startheston = time.time()
    process = Popen(["../bin/heston-pricing", "--iters", str(100*i)], stdout=PIPE)
    (_, _) = process.communicate(); _ = process.wait()
    stopheston = time.time()
    

    # Asian BNS
    # BNSAsianCall(s0, r, T, K, rho, k, theta, dzeta, alpha, c, lambda, phi)
    # BNSAsianCall asset(50, 0.05, 1, 44, 0.0, 1, 0.00, 0, 0.5, 0.01, 1, -1.0); 
    startbns = time.time()
    process = Popen(["../bin/heston-pricing", "--model", "bns",  "--meanrevert", "0", "--correlation", "0", "--elasticity", "1",  "--iters", str(100*i)], stdout=PIPE)
    (_, _) = process.communicate(); _ = process.wait()
    stopbns = time.time()
    x.append(100*i)
    yheston.append(stopheston - startheston)
    ybns.append(stopbns - startbns)

plt.figure()
plt.title("Comparison of the computation times")
plt.plot(x, yheston, label="Asian call in Heston model")
plt.plot(x, ybns, label="Asian call in BNS model")
plt.xlabel("Number of iterations")
plt.ylabel("Time to compute (seconds)")
plt.legend()
plt.show()