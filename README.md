# Estimating Hofer-Zehnder Capacities

This repository contains development of earlier work by [1] and the Symplectic Capacities group from the Computational Symplectic Topology workshop at Tel-Aviv University and ICERM [2], of which I was a participant.  Work by the Symplectic Capacities group is available in the github repo [3].

The contents of the repository are currently as follows:

## HZCapacityEstimator.py 

Defines a class that estimates the Hofer-Zehnder capacity of a unit-ball in R^2n. The algorithm is a vectorized version of the constrained gradient-descent algorithm provided in [1].  In the algorithm, loops are approximated by piecewise linear paths with m linear segments.

Usage:

- After importing '''HZCapacityEstimator''', initialize an instance of the class by calling

''' estimator = HZCapacityEstimator(n,m)'''

where n and m are as described above (n is half the desired dimension and m is the number of segments in the PL paths).  Note that m must be at least 3. For a good approximation, we recommend setting m to be approximately 1000.

- Perform gradient descent by calling the function 

''' estimator.estimate(iterations = 100, epsilon = 10.**(-12))'''

iterations is the maximum number of gradient-descent steps you wish to take, and epsilon is the tolerance for early stopping.


##




## References

[1] G\"oing-Jaeschke. Numerical Analysis in Hamiltonian Dynamics. PhD Thesis.

[2] http://www.math.tau.ac.il/~ostrover/Workshop/CST/Workshop-CST.html

[3] https://github.com/itamarr/symp-cap