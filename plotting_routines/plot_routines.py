import numpy as np
import matplotlib.pyplot as plt

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
    import booz_xform as bx
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
    plt.close()
    import  matplotlib       as      mpl
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



def plot_3D(vmec):
    plt.close()
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    rmnc    = vmec.wout.rmnc.T
    zmns    = vmec.wout.zmns.T
    bmnc    = vmec.wout.bmnc.T
    ns      = vmec.wout.ns
    xn = vmec.wout.xn
    xm = vmec.wout.xm
    xn_nyq = vmec.wout.xn_nyq
    xm_nyq = vmec.wout.xm_nyq
    rmns = 0*rmnc
    zmnc = 0*rmnc
    bmns = 0*bmnc
    nmodes = len(xn)
    fig = plt.figure("3D Surface Plot")
    ntheta = 100
    nzeta = int(200)
    theta1D = np.linspace(0,2*np.pi,num=ntheta)
    zeta1D = np.linspace(0,2*np.pi*3/4,num=nzeta)
    zeta2D, theta2D = np.meshgrid(zeta1D,theta1D)
    iradius = ns-1
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

    fig.patch.set_facecolor('white')
    ax = fig.add_subplot(111, projection='3d')
    #ax = fig.gca(projection='3d',azim=0, elev = 45/2)
    ax._axis3don = False
    ax.plot_surface(X, Y, Z, facecolors = cm.turbo(B_rescaled), rstride=1, cstride=1, antialiased=False)
    ax.auto_scale_xyz([X.min(), X.max()], [X.min(), X.max()], [X.min(), X.max()])
    plt.show()



def plot_cross_section(vmec):
    import math
    # import a lot of stuff 
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
    iradius = ns-1
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

    
    fig = plt.figure("Poincare Plots",figsize=(14,7))
    fig.patch.set_facecolor('white')

    numCols = 3
    numRows = 2
    plotNum = 1

    plt.subplot(numRows,numCols,plotNum)
    #plt.subplot(1,1,1)
    plotNum += 1
    for ind in range(nzeta):
        plt.plot(R[:,ind], Z[:,ind], '-')
        plt.plot(R[:,ind], Z[:,ind], '-')
        plt.plot(R[:,ind], Z[:,ind], '-')
        plt.plot(R[:,ind], Z[:,ind], '-')
    plt.gca().set_aspect('equal',adjustable='box')
    #plt.legend(fontsize='x-small')
    plt.xlabel('R')
    plt.ylabel('Z')


    ntheta = 200
    nzeta = 5
    nradius = 10
    theta = np.linspace(0,2*np.pi,num=ntheta)
    zeta = np.linspace(0,2*np.pi/nfp/2,num=nzeta,endpoint=True)
    iradii = np.linspace(0,ns-1,num=nradius).round()
    iradii = [int(i) for i in iradii]
    R = np.zeros((ntheta,nzeta,nradius))
    Z = np.zeros((ntheta,nzeta,nradius))
    for itheta in range(ntheta):
        for izeta in range(nzeta):
            for iradius in range(nradius):
                for imode in range(nmodes):
                    angle = xm[imode]*theta[itheta] - xn[imode]*zeta[izeta]
                    R[itheta,izeta,iradius] = R[itheta,izeta,iradius] + rmnc[iradii[iradius],imode]*math.cos(angle) \
                                                                    + rmns[iradii[iradius],imode]*math.sin(angle)
                    Z[itheta,izeta,iradius] = Z[itheta,izeta,iradius] + zmns[iradii[iradius],imode]*math.sin(angle) \
                                                                    + zmnc[iradii[iradius],imode]*math.cos(angle)
    for izeta in range(nzeta):
        plt.subplot(numRows,numCols,plotNum)
        plotNum += 1
        for iradius in range(nradius):
            plt.plot(R[:,izeta,iradius], Z[:,izeta,iradius], '-')
        plt.plot(Raxis[izeta],Zaxis[izeta],'xr')
        plt.gca().set_aspect('equal',adjustable='box')
        plt.xlabel('R')
        plt.ylabel('Z')
        plt.title('zeta = '+str(zeta[izeta]))

    #maximizeWindow()

    plt.tight_layout()
    plt.show()



# plot info 
def plot_vmec_data(vmec):
    import math
    # import a lot of stuff 
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
    iradius = ns-1
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
        
    


    xLabel = r'$s = \psi_N$'
    fig = plt.figure("VMEC Data",figsize=(14,7))
    fig.patch.set_facecolor('white')

    numCols = 3
    numRows = 3
    plotNum = 1

    plt.subplot(numRows,numCols,plotNum)
    plotNum += 1
    #plt.plot(phi, iotaf, '.-',label='iotaf')
    #plt.plot(phi_half, iotas[1:],'.-',label='iotas')
    plt.plot(s, iotaf, '.-',label='iotaf')
    plt.plot(s_half, iotas[1:],'.-',label='iotas')
    plt.legend(fontsize='x-small')
    plt.xlabel(xLabel)

    plt.subplot(numRows,numCols,plotNum)
    plotNum += 1
    #plt.plot(phi, presf, '.-',label='presf')
    #plt.plot(phi_half, pres[1:], '.-',label='pres')
    plt.plot(s, presf, '.-',label='presf')
    plt.plot(s_half, pres[1:], '.-',label='pres')
    plt.legend(fontsize='x-small')
    plt.xlabel(xLabel)

    plt.subplot(numRows,numCols,plotNum)
    plotNum += 1
    #plt.plot(phi_half, buco[1:], '.-',label='buco')
    plt.plot(s_half, buco[1:], '.-',label='buco')
    plt.title('buco')
    plt.xlabel(xLabel)

    plt.subplot(numRows,numCols,plotNum)
    plotNum += 1
    #plt.plot(phi_half, bvco[1:], '.-',label='bvco')
    plt.plot(s_half, bvco[1:], '.-',label='bvco')
    plt.title('bvco')
    plt.xlabel(xLabel)

    plt.subplot(numRows,numCols,plotNum)
    plotNum += 1
    #plt.plot(phi, jcuru, '.-',label='jcuru')
    plt.plot(s, jcuru, '.-',label='jcuru')
    plt.title('jcuru')
    plt.xlabel(xLabel)

    plt.subplot(numRows,numCols,plotNum)
    plotNum += 1
    #plt.plot(phi, jcurv, '.-',label='jcurv')
    plt.plot(s, jcurv, '.-',label='jcurv')
    plt.title('jcurv')
    plt.xlabel(xLabel)

    plt.subplot(numRows,numCols,plotNum)
    plotNum += 1
    # if 'power_series' in pcurr_type:
    # 	ac_profile = phi*0.0
    # 	for i in range(len(ac)):
    # 		ac_profile += ac[i]*(s**i)
    # 	plt.plot(s, ac_profile, '.-')
    # else:
    # 	mask = (ac_aux_s >= 0)
    # 	plt.plot(ac_aux_s[mask], ac_aux_f[mask],'.-')
    # plt.title('ac profile')
    # plt.xlabel(xLabel)
    s_plot_ignore=0.3
    plt.plot(s[int(s_plot_ignore*len(s)):-2],DMerc[int(s_plot_ignore*len(s)):-2],'.-')
    plt.title('DMerc')
    plt.xlabel(xLabel)
    plt.tight_layout()
    plt.show()