import numpy as np
from simsopt.util.mpi                   import MpiPartition
from simsopt.mhd                        import Vmec, vmec_fieldlines
from mpi4py                             import MPI
from plot_routines                      import *

mpi = MpiPartition()
mpi.write()

# please choose a radial coordinate (s = (r/a)**2), and a field-line label (alpha)
s_val       = 1.0
alpha_val   = 0.0


# let us first generate the vmec files from the input files
vmec_QA     = Vmec(filename="wout_precise_QA.nc",mpi=mpi,verbose=True)

# let us also generate some magnetic field lines
# first generate a grid of poloidal points on which to generate the fieldline
n_pol       = 2.0
theta_grid  = np.linspace(-n_pol*np.pi,n_pol*np.pi,1000)
fl_QA       = vmec_fieldlines(vmec_QA, s_val, alpha_val, theta1d=theta_grid)

# plot the flux-surface and field line on it
plot_surface_and_fl(vmec_QA,fl_QA,s_val,transparant=False,trans_val=0.9,title='')


