# Estimating Hofer-Zehnder Capacities

This repository contains development of earlier work by [1] and the Symplectic Capacities group from the Computational Symplectic Topology workshop at Tel-Aviv University and ICERM [[2]](http://www.math.tau.ac.il/~ostrover/Workshop/CST/Workshop-CST.html), of which I was a participant.  Work by the Symplectic Capacities group is available in MATLAB code in the github repo [[3]](https://github.com/itamarr/symp-cap).

Hofer-Zehder capacities of convex sets in R^2n are an example of an invariant of symplectic manifolds (another example is Gromov width). These invariants are generally very difficult to compute or even estimate, and there are many open conjectures, even about the Hofer-Zehder capacity.

The contents of the repository are currently as follows:

## HZCapacityEstimator.py 

Defines a class that estimates the Hofer-Zehnder capacity of a strictly convex open subset of R^2n. The algorithm is a vectorized numpy implementation of the constrained gradient-descent algorithm provided in [1].  The constraints are non-linear and the algorithm uses a version of projected gradients that solves the local KKT conditions.  In the algorithm, smooth loops in R^2n are approximated by piecewise linear paths with m linear segments.

Usage of ```HZCapacityEstimator.py``` is demonstrated in the notebook ```Examples.ipynb```

### Initializing ```HZCapacityEstimator```

After importing ```HZCapacityEstimator```, initialize an instance of the class by calling 

```estimator = HZCapacityEstimator(n,m,H)```

#### Arguments

- ```n``` is half the desired dimension.

- ```m``` is the number of segments in the PL paths).  Note that m must be at least 3. For a good approximation, we recommend setting m to be approximately 500-1000.

- ```H``` is a function 

### Running the ```HZCapacityEstimator.estimate```

Once you have defined an estimator object, as above, you can run the gradient-descent algorithm by calling 

```estimator.estimate(iterations = 100, epsilon = 10.**(-12), VOCAL = True)```

#### Arguments

- ```iterations``` this optional argument sets the maximum number of gradient descent steps.

- ```epsilon``` this optional argument is the tolerance for early stopping.

- ```VOCAL``` this optional argument can be used to mute log messages that the algorithm prints to screen. 

#### Returns

The function returns the estimated HZ- capacity.


## BallHZCapacityEstimator.py 

Same as ```HZCapacityEstimator.py``` but hardcoded with Hamiltonian H defining the unit ball in R^2n. It is known that the HZ-capacity of a unit ball is pi.

## To do 

- Expand the algorithm in ```TropicalHZCapacityEstimator.py``` to estimate capacities of convex polytopes. 
- Add visualization code.


## References

[1] G\"oing-Jaeschke. Numerical Analysis in Hamiltonian Dynamics. PhD Thesis, ETH Zurich.

[2] http://www.math.tau.ac.il/~ostrover/Workshop/CST/Workshop-CST.html

[3] https://github.com/itamarr/symp-cap