import numpy as np
from simsopt.util.mpi                   import MpiPartition
from simsopt.mhd                        import Vmec, vmec_fieldlines
from mpi4py                             import MPI
import matplotlib.pyplot as plt

mpi = MpiPartition()
mpi.write()


### In this script we simply run a vmec input file to get wout files.
### Wout files can again be read by vmec and used to generate fieldlines and various 
### other quantities.


# give input file name
filename    = "input.precise_QA"
# run vmec. verbose command shows progress of the VMEC run.
vmec        = Vmec(filename=filename,mpi=mpi,verbose=True)
