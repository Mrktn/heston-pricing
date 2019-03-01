#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
    print("No file supplied")
    sys.exit()

with open(sys.argv[1]) as f:
    l = []
    for x in f:
        l.append(float(x))

    plt.hist(l)
    plt.title("Average = " + str(round(np.mean(l), 3)) + ", std = " + str(round(np.std(l), 3)))
    plt.show()
