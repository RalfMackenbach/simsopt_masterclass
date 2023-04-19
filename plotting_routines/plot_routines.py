import numpy as np


def plot_surface_and_fl(vmec,fl,s_val,transparant=False,trans_val=0.9,title=''):
    import  math
    from    matplotlib          import  cm
    from    mayavi import mlab
    phi = vmec.wout.phi
    iotaf = vmec.wout.iotaf
    presf = vmec.wout.presf
    iotas = vmec.wout.iotas
    pres = vmec.wout.pres
    ns = vmec.wout.ns
    nfp = vmec.wout.nfp
    xn = vmec.wout.xn
    xm = vmec.wout.xm
    xn_nyq = vmec.wout.xn_nyq
    xm_nyq = vmec.wout.xm_nyq
    rmnc = vmec.wout.rmnc.T
    zmns = vmec.wout.zmns.T
    bmnc = vmec.wout.bmnc.T
    raxis_cc = vmec.wout.raxis_cc
    zaxis_cs = vmec.wout.zaxis_cs
    buco = vmec.wout.buco
    bvco = vmec.wout.bvco
    jcuru = vmec.wout.jcuru
    jcurv = vmec.wout.jcurv
    lasym = vmec.wout.lasym
    iradius = int(np.floor(s_val*(ns-1)))

    ac_aux_s = vmec.wout.ac_aux_s
    ac_aux_f = vmec.wout.ac_aux_f

    mpol = vmec.wout.mpol
    ntor = vmec.wout.ntor
    Aminor_p = vmec.wout.Aminor_p
    Rmajor_p = vmec.wout.Rmajor_p
    aspect = vmec.wout.aspect
    betatotal = vmec.wout.betatotal
    betapol = vmec.wout.betapol
    betator = vmec.wout.betator
    betaxis = vmec.wout.betaxis
    ctor = vmec.wout.ctor
    DMerc = vmec.wout.DMerc
    gmnc = vmec.wout.gmnc

    if lasym == 1:
        rmns = vmec.wout.rmns
        zmnc = vmec.wout.zmnc
        bmns = vmec.wout.bmns
        raxis_cs = vmec.wout.raxis_cs
        zaxis_cc = vmec.wout.zaxis_cc
    else:
        rmns = 0*rmnc
        zmnc = 0*rmnc
        bmns = 0*bmnc
        raxis_cs = 0*raxis_cc
        zaxis_cc = 0*raxis_cc
    try:
        ac = vmec.wout.ac
    except:
        ac = []
    try:
        pcurr_type = vmec.wout.pcurr_type
    except:
        pcurr_type = ""
    nmodes = len(xn)
    s = np.linspace(0,1,ns)
    s_half = [(i-0.5)/(ns-1) for i in range(1,ns)]
    phiedge = phi[-1]
    phi_half = [(i-0.5)*phiedge/(ns-1) for i in range(1,ns)]
    ntheta = 200
    nzeta = 8
    theta = np.linspace(0,2*np.pi,num=ntheta)
    zeta = np.linspace(0,2*np.pi/nfp,num=nzeta,endpoint=False)
    R = np.zeros((ntheta,nzeta))
    Z = np.zeros((ntheta,nzeta))
    for itheta in range(ntheta):
        for izeta in range(nzeta):
            for imode in range(nmodes):
                angle = xm[imode]*theta[itheta] - xn[imode]*zeta[izeta]
                R[itheta,izeta] = R[itheta,izeta] + rmnc[iradius,imode]*math.cos(angle) + rmns[iradius,imode]*math.sin(angle)
                Z[itheta,izeta] = Z[itheta,izeta] + zmns[iradius,imode]*math.sin(angle) + zmnc[iradius,imode]*math.cos(angle)

    Raxis = np.zeros(nzeta)
    Zaxis = np.zeros(nzeta)
    for izeta in range(nzeta):
        for n in range(ntor+1):
            angle = -n*nfp*zeta[izeta]
            Raxis[izeta] += raxis_cc[n]*math.cos(angle) + raxis_cs[n]*math.sin(angle)
            Zaxis[izeta] += zaxis_cs[n]*math.sin(angle) + zaxis_cc[n]*math.cos(angle)
    ntheta = 100
    nzeta = int(200)
    theta1D = np.linspace(0,2*np.pi,num=ntheta)
    zeta1D = np.linspace(0,2*np.pi,num=nzeta)
    zeta2D, theta2D = np.meshgrid(zeta1D,theta1D)
    R = np.zeros((ntheta,nzeta))
    Z = np.zeros((ntheta,nzeta))
    B = np.zeros((ntheta,nzeta))
    for imode in range(nmodes):
        angle = xm[imode]*theta2D - xn[imode]*zeta2D
        R = R + rmnc[iradius,imode]*np.cos(angle) + rmns[iradius,imode]*np.sin(angle)
        Z = Z + zmns[iradius,imode]*np.sin(angle) + zmnc[iradius,imode]*np.cos(angle)

    for imode in range(len(xn_nyq)):
        angle = xm_nyq[imode]*theta2D - xn_nyq[imode]*zeta2D
        B = B + bmnc[iradius,imode]*np.cos(angle) + bmns[iradius,imode]*np.sin(angle)

    X = R * np.cos(zeta2D)
    Y = R * np.sin(zeta2D)
    # Rescale to lie in [0,1]:
    B_rescaled = (B - B.min()) / (B.max() - B.min())


    X_coord = (fl.R).flatten() * (fl.cosphi).flatten()
    Y_coord = (fl.R).flatten() * (fl.sinphi).flatten()
    Z_coord = (fl.Z).flatten()

    from mayavi import mlab
    fig = mlab.figure(size=(1000, 1000))

    colors = np.asarray(cm.binary(B_rescaled))[:,:,1]
    surf3 = mlab.mesh(X, Y, Z, scalars=B,colormap='coolwarm',vmin=np.amin(B),vmax=np.amax(B))
    line = mlab.plot3d(X_coord, Y_coord, Z_coord,color=(0.0,0.0,0.0),tube_radius=0.03) #(0, 0.7, 0.23)
    gridpoints = len(X_coord)
    point = mlab.points3d(X_coord[int(gridpoints/2)], Y_coord[int(gridpoints/2)], Z_coord[int(gridpoints/2)], color=(1.0,1.0,1.0),scale_factor=0.2)

    if transparant == True:
        surf3.actor.property.opacity = trans_val

    if title!='':
       mlab.savefig('3D_plot_'+title+'.png', figure=fig)



    mlab.show()
    mlab.close(all=True)



def plot_boozer(wout_name,s_val):
    import matplotlib.pyplot as plt
    import booz_xform as bx
    print(s_val)
    b1 = bx.Booz_xform()
    b1.read_wout(wout_name)
    b1.compute_surfs = [s_val]
    b1.mboz = 40
    b1.nboz = 40
    b1.run()
    #fig = plt.figure()
    bx.surfplot(b1, js=0,  fill=False, ncontours=40, ntheta=200,nphi=200)
    plt.show()


def plot_geom(fl):
    import  matplotlib       as      mpl
    import matplotlib.pyplot as plt
    font = {'family': 'sans-serif',
        'weight': 'normal',
        'size': 10}

    mpl.rc('font', **font)

    theta = (fl.theta_pest).flatten()/np.pi

    fig, ax = plt.subplots(2, 4, figsize=(10, 5),tight_layout=True)
    # plot magnetic field strength
    ax[0,0].plot(theta,(fl.modB).flatten())
    ax[0,0].set_xlabel(r'$\theta/\pi$')
    ax[0,0].set_ylabel(r'$|B|$')
    # plot nabla psi . nabla psi
    ax[0,1].plot(theta,(fl.grad_psi_dot_grad_psi).flatten())
    ax[0,1].set_xlabel(r'$\theta/\pi$')
    ax[0,1].set_ylabel(r'$\nabla \psi \cdot \nabla \psi$')
    # plot nabla alpha . nabla alpha
    ax[0,2].plot(theta,(fl.grad_alpha_dot_grad_alpha).flatten())
    ax[0,2].set_xlabel(r'$\theta/\pi$')
    ax[0,2].set_ylabel(r'$\nabla \alpha \cdot \nabla \alpha$')
    # plot d l / d theta
    dtdl = (fl.B_sup_theta_pest/fl.modB).flatten()
    ax[0,3].plot(theta,np.abs(1/dtdl))
    ax[0,3].set_xlabel(r'$\theta/\pi$')
    ax[0,3].set_ylabel(r'$\mathrm{d \ell} / \mathrm{d} \theta$')
    # plot radial projection of gradient drift 
    ax[1,0].plot(theta,(fl.B_cross_grad_B_dot_grad_psi).flatten())
    ax[1,0].set_xlabel(r'$\theta/\pi$')
    ax[1,0].set_ylabel(r'$\mathbf{B} \times \nabla B \cdot \nabla \psi $')
    # plot radial projection of curvature drift 
    ax[1,1].plot(theta,(fl.B_cross_kappa_dot_grad_psi).flatten())
    ax[1,1].set_xlabel(r'$\theta/\pi$')
    ax[1,1].set_ylabel(r'$\mathbf{B} \times \mathrm{\kappa} \cdot \nabla \psi $')
    # plot binormal projection of gradient drift 
    ax[1,2].plot(theta,(fl.B_cross_grad_B_dot_grad_alpha).flatten())
    ax[1,2].set_xlabel(r'$\theta/\pi$')
    ax[1,2].set_ylabel(r'$\mathbf{B} \times \nabla B \cdot \nabla \alpha $')
    # plot radial projection of curvature drift 
    ax[1,3].plot(theta,(fl.B_cross_kappa_dot_grad_alpha).flatten())
    ax[1,3].set_xlabel(r'$\theta/\pi$')
    ax[1,3].set_ylabel(r'$\mathbf{B} \times \mathrm{\kappa} \cdot \nabla \alpha $')

    plt.show()