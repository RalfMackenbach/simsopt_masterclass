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
    iota = vmec.iota_axis()
    print('iota on axis = ', iota)
    pen = iota-t
    return pen


# quasisymmetry target
def qs(vmec,t=0.0,M=1,N=0):
    """
    A quasi-symmetry target. Is zero in exactly QS systems.
    M and N is helicty. M=1, N=0 is the quasi-axisymmetry.
    """
    qs = QuasisymmetryRatioResidual(vmec,
                                    np.arange(0, 1.01, 0.1),  # Radii to target
                                    helicity_m=M, helicity_n=N) 
    pen = np.linalg.norm(qs.residuals() - t)
    print("Quasi-symmetry penalty =", pen)
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
    s_q = (fl.shat).flatten()[0]

    pen = s_q - t
    
    ##########################################################################
    ##########################################################################
    ##########################################################################
    print("My target penalty =", pen)
    return pen