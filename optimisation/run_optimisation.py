
#!/usr/bin/env python

import os
import numpy as np
from simsopt.util import MpiPartition
from simsopt.mhd import Vmec
from simsopt.objectives import LeastSquaresProblem
from simsopt.solve import least_squares_mpi_solve
from simsopt import make_optimizable
from target_functions import *

"""
This script run the optimisation using various targets defined
in target_functions.py. If you want a different stellarator 
periodicity (e.g. Wendelstein 7-X has NFP=5), you need to change
NFP in input.tok to your chosen value.
"""

print("Running optimisation...")
print("=======================")



mpi = MpiPartition()
mpi.write()
# Create initial condition for the optimisation.
# We will start with a circular tokamak without current.
vmec = Vmec(filename="input.tok",mpi=mpi,verbose=False)
surf = vmec.boundary

# make list in which we will store the targets
make_tuples = []


# Define and append targets to make_tuples. The two values after the target are the wanted value, and the weight of the target.

# aspect_target
make_tuples.append((make_optimizable(aspect_target, vmec,8.0,True).J, 0.0, 1.0))

# iota target
make_tuples.append((make_optimizable(iota0_target, vmec, 0.5).J, 0.0, 1.0))

# quasi-symmetry target
make_tuples.append((make_optimizable(qs, vmec,0.0).J, 0.0, 1.0))

# my target
make_tuples.append((make_optimizable(my_fl_target, vmec, 0.5).J, 0.0, 1.0))


# Define objective function
prob = LeastSquaresProblem.from_tuples(make_tuples)


# simsopt really excels in dynamic optimization,
# where the resolution of the problem changes
# as we are iterating. We solve the problem
# for increasingly higher resolution (mpol, ntor).
for step in range(3):
    max_mode = step + 1

    # VMEC's mpol & ntor will be 2, 3, 4, 5:
    vmec.indata.mpol = 2 + step
    vmec.indata.ntor = vmec.indata.mpol

    if mpi.proc0_world:
        print("Beginning optimization with max_mode =", max_mode, \
              ", vmec mpol=ntor=", vmec.indata.mpol, \
              ". Previous vmec iteration = ", vmec.iter)

    # Define parameter space:
    surf.fix_all()
    surf.fixed_range(mmin=0, mmax=max_mode, 
                     nmin=-max_mode, nmax=max_mode, fixed=False)
    surf.fix("rc(0,0)")  # Major radius

    # For the test to run quickly, we stop after the first function
    # evaluation, by passing max_nfev=1 to scipy.optimize. For a
    # "real" optimization, remove the max_nfev parameter below.
    least_squares_mpi_solve(prob, mpi, grad=True)

    # Preserve the output file from the last iteration, so it is not
    # deleted when vmec runs again:
    vmec.files_to_delete = []

    if mpi.proc0_world:
        print(f"Done optimization with max_mode ={max_mode}. "
              f"Final vmec iteration = {vmec.iter}")

print("Good bye")