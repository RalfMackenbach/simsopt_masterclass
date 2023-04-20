import numpy as np
from simsopt.util.mpi                   import MpiPartition
from simsopt.mhd                        import Vmec, vmec_fieldlines
from mpi4py                             import MPI
from plot_routines                      import *
import matplotlib.pyplot as plt

mpi = MpiPartition()
mpi.write()

# please choose a radial coordinate (s = (r/a)**2), and a field-line label (alpha)
s_val       = 0.9
alpha_val   = 0.0


# read in your wout file here.
filename    = "wout_precise_QA.nc"
vmec     = Vmec(filename=filename,mpi=mpi,verbose=True)

# let us also generate some magnetic field lines
# first generate a grid of poloidal points on which to generate the fieldline
n_pol       = 1.0
theta_grid  = np.linspace(-n_pol*np.pi,n_pol*np.pi,1000)
fl       = vmec_fieldlines(vmec, s_val, alpha_val, theta1d=theta_grid)

# plot the flux-surface and field line on it
plot_surface_and_fl(vmec,fl,s_val,transparant=False,trans_val=0.9,title='')

# plot surface alone
plot_3D(vmec)

# plot the Boozer surface
# plot_boozer takes as input the filename and the radial surface NUMBER.
# Let us first find this number:
ns = len(vmec.s_full_grid)
s_num = int(s_val*ns)
plot_boozer(filename,s_num)

# plot geometric quantities of the field line
plot_geom(fl)

# plot cross sections
plot_cross_section(vmec)

# plot data
plot_vmec_data(vmec)