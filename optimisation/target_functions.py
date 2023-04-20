import  numpy                       as      np
from simsopt.mhd.vmec               import  Vmec
from simsopt.mhd.vmec_diagnostics   import  vmec_fieldlines
from simsopt.mhd                    import  QuasisymmetryRatioResidual

# volume target
def volume_target(vmec,vol,one_sided=True):
    """
    Targets the volume.
    If one sided is False, the penalty is only zero iff the volume is exactly equal to vol.
    If one sided is True, the penalty is zero iff the volume is greater than vol.
    """
    vmec.run()
    vmec_volume = vmec.volume()
    targ_volume = vol
    print("V_vmec/V_want = ",vmec_volume/targ_volume)
    pen = 1-vmec_volume/targ_volume
    if one_sided == True:
        pen = np.max([0,pen])
    return pen


# aspect ratio target
def aspect_target(vmec,t=2.0,one_sided=True):
    """
    Aspect ratio target.
    If one sided is False, the penalty is only zero iff the aspect ratio is exactly equal to t.
    If one sided is True, the penalty is zero iff the aspect ratio is greater than t.
    """
    vmec.run()
    aspect_ratio = vmec.aspect()
    print('Aspect ratio = ', aspect_ratio)
    pen = aspect_ratio-t
    if one_sided == True:
        pen = np.max([0,pen])
    return pen


# rotational transform on axis target
def iota0_target(vmec,t=0.5):
    """
    Rotational transform on axis target.
    """
    vmec.run()
    iota = np.abs(vmec.iota_axis())
    print('iota on axis = ', iota)
    pen = iota-t
    return pen


# quasisymmetry target
def qs(vmec,t=0.0,M=1,N=1):
    """
    A quasi-symmetry target. Is zero in exactly QS systems.
    M and N is helicty. M=1, N=0 is the quasi-axisymmetry,
    and M=1 N=1 is the quasi-helical symmetry.
    """
    qs = QuasisymmetryRatioResidual(vmec,
                                    np.arange(0, 1.01, 0.1),  # Radii to target
                                    helicity_m=M, helicity_n=N) 
    pen = np.sum(np.abs(qs.residuals() - t))
    print("Quasi-symmetry penalty =", pen)
    return pen


# penalize mirror ratio, i.e. makes sure the variation in magnetic field is not too large
def mirror_ratio_target(vmec,t=0.21):
    """
    Adopted from A. Goodman. Calculates mirror near axis. Penalizes iff the target value is exceeded.
    vmec        -   the vmec object
    t           -   target value
    """
    vmec.run()
    xm_nyq = vmec.wout.xm_nyq
    xn_nyq = vmec.wout.xn_nyq
    bmnc = vmec.wout.bmnc.T
    bmns = 0*bmnc
    nfp = vmec.wout.nfp

    Ntheta = 100
    Nphi = 100
    thetas = np.linspace(0,2*np.pi,Ntheta)
    phis = np.linspace(0,2*np.pi/nfp,Nphi)
    phis2D,thetas2D=np.meshgrid(phis,thetas)
    b = np.zeros([Ntheta,Nphi])
    for imode in range(len(xn_nyq)):
        angles = xm_nyq[imode]*thetas2D - xn_nyq[imode]*phis2D
        b += bmnc[1,imode]*np.cos(angles) + bmns[1,imode]*np.sin(angles)
    Bmax = np.max(b)
    Bmin = np.min(b)
    m = (Bmax-Bmin)/(Bmax+Bmin)
    print("Mirror = ",m)
    pen = np.max([0,1-t/m])
    return pen


# my target function
def my_fl_target(vmec,t,s_val=0.5,alpha_val=0.0,n_pol=1.0,theta_res=1000):
    theta_grid  = np.linspace(-n_pol*np.pi,n_pol*np.pi,theta_res)
    vmec.run()
    fl          = vmec_fieldlines(vmec, s_val, alpha_val, theta1d=theta_grid)

    ##########################################################################
    ### construct your target using the magnetic field line object fl here ###
    ##########################################################################
    
    # as an example, here's how to target the magnetic shear at s=0.5
    # shear is ill-defined at q=0, so let do r*dq/dr instead.
    s_q = np.abs(((fl.iota)*(fl.shat)).flatten())
    
    pen = s_q - t

    ##########################################################################
    ##########################################################################
    ##########################################################################
    print("My target penalty =", pen)
    return pen


### Feel free to add more target functions!