#!/usr/bin/python
from subprocess import Popen, PIPE
import numpy as np
import matplotlib.pyplot as plt 
plt.style.use('seaborn-whitegrid')

strikes = range(30, 81, 2)
vals = np.array([])
stds = np.array([])

N = 100

for K in strikes:
    print(K)
    l = np.array([])
    lb = np.array([])
    for i in range(N):
        # Asian call heston
        #process = Popen(["../bin/heston-pricing", "--strike", str(K), "-t", "asian"], stdout=PIPE)
        # Asian call BNS with P&P parameters
        process = Popen(["../bin/heston-pricing", "--model", "bns",  "--meanrevert", "0", "--correlation", "0", "--elasticity", "1", "--strike", str(K), "--volofvar", "0"], stdout=PIPE)
        (output, _) = process.communicate(); _ = process.wait()
        l = np.append(l, float(output))
    vals  = np.append(vals, np.mean(l))
    stds  = np.append(stds, np.std(l))

plt.close('all')
plt.plot(strikes, vals, alpha=0.7, color='black', label="Price")
plt.fill_between(strikes, vals - stds, vals + stds, color='black', alpha=0.2, label="1 standard deviation")

plt.title("Asian call option price dependency to the strike (BNS model)")

plt.xlabel("$K$")
plt.ylabel("$C(K)$")
plt.legend()
plt.show()

print(list(strikes))
print(list(vals))
print(list(stds))