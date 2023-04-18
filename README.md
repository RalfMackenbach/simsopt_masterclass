# simsopt_masterclass
A repo for the simsopt part of the stellarator masterclass at TU/e

Install simsopt via the installation instructions on the website: https://simsopt.readthedocs.io/en/latest/installation.html .
To keep your device tidy, I would recommend using a virtual environment manager such as anaconda. If you have an M1 chip, please follow: https://github.com/hiddenSymmetries/simsopt/wiki/Mac-M1-installation

Please also install the "optional" packages: mpi4py, booz_xform, and VMEC.
mpi4py can readily be installed from your package manager (e.g. ```conda install mpi4py```).
booz_xform can simply be installed with ```pip install booz_xform``` within your virtual environment (in some cases netCDF errors may occur - please do ```conda install netcdf-fortran``` in such cases)

VMEC can be more challenging to get up and running. I recommend following the GitHub page: https://github.com/hiddenSymmetries/VMEC2000 . Important, if you have a failed installation remove the folder ```_skbuild``` in the VMEC2000 repo before trying to rebuild! If you have all installed, you are ready to play around with ```simsopt```!

For some plotting routines I'd recommend installing mayavi, e.g. ```conda install mayavi```
