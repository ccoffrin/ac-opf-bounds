# Bound Tightening for AC Optimal Power Flow

This code provides a minimalist python+gurobi based implementation for computing tight bounds on AC-OPF test cases in the [matpower](http://www.pserc.cornell.edu/matpower/) format.  Formally, the algorithm computes the "minimal continuous constraint relaxation network" of the AC-OPF problem.  Namely, the finds new bounds on the bus voltage magnitudes and the line phase angle differences (i.e. pad) that do not remove any feasible solution to AC-OPF constraints.  The method is independent of the objective function under consideration and is applicable with a wide range of power system applications using the operational constraint from the AC-OPF problem.  

Details on the theory and correctness of the algorithm can be found in this [paper](http://link.springer.com/chapter/10.1007%2F978-3-319-23219-5_4).  The method is based on the [Quadratic Convex (QC) relaxation](http://www.optimization-online.org/DB_HTML/2013/09/4057.html) of AC power flows, and includes the "lifted nonlinear cuts" from this [report](http://arxiv.org/abs/1512.04644).


## Installation

The primary entry point is the "compute-bounds.py" scipt.  It can be run in place without installation.  However, it does require gurobi to be installed so that it can imported as a python library, i.e. "import gurobipy".  See [gurobi](http://gurobi.com) for detailed installation instructions.


## Usage

The standard way to use this script is to provide a matpower case file as follows,
```
compute-bounds.py case.m
```
This will solve _2(n+l)_ QC relaxations in a number of rounds, where _n_ is the number nodes in the network and _l_ is the number of lines.  Note that each round may take several minutes for cases with over 100 buses.  The progress of each round is summarized by an "ITER_DATA" line.  Once a fixpoint is reached (i.e. no further tightening can be made) the script completes by printing out tables of the updated voltage and angle bounds as well as a summary data in a "SMRY_DATA" line.

The default computation can be modified using command line arguments.  For example, "--large" can be used for increasing numerical stability on cases with over 1000 buses and "-output" can be used to see more details of the computation.  Use "--help" to see a complete list of options.


## Comparison to Other Methods

The [paper](http://link.springer.com/chapter/10.1007%2F978-3-319-23219-5_4) that proposed this algorithm observed that bound tightening was a powerful pre-processing step for solving AC-OPF problems.  Notably, this bound tightening procedure can make the Convex Quadratic (QC) relaxation stronger than the well-known Semi-Definite Relaxation (SDP).  This is appealing becouse quadratic programming solvers, such as gurobi, are more scalable and reliable than SDP solvers. 

Running this script can consume a significant amount of time for cases with over 100 buses.  However, it is important to note that this is a single threaded implementation of a highly parallelizable computation.  The "parallel time" field in the summary output suggests the runtime of the algorithm with maximum parallelization.  The ideal parallel runtime is less than 5 minutes on all publicly available test cases.


## Development

Community-driven development and enhancement of this script are welcome and encouraged.  Please fork this repository and share your contributions to the master with pull requests.


## Acknowledgments

Although this code was developed by Carleton Coffrin, the theory behind this algorithm was done in collaboration with Hassan L. Hijazi, and Pascal Van Hentenryck.


## License

MIT, see LICENSE file for the details.
