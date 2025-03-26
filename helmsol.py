from ngsolve import *
import numpy as np
ngsglobals.msg_level = 0

def farf2d_penetrable(u,mesh,farp,kappa,porder):
    # Input
    # u - the near field solution computed by helmsol
    # mesh - the NGSpy mesh
    # farp - parameters for the far field pattern
    # kappa - the wavenumber
    # porder - the order of the polynmials for the FE calculation
    #
    # Output
    # uinf - the far field pattern
    # phi - angles for directions where uinf is calculated
    nphi=farp["n"]
    cphi=farp["cent"]
    appphi=farp["app"]
    nv=specialcf.normal(mesh.dim)
    phi=np.zeros(nphi)
    uinf=np.zeros(nphi,dtype=complex)
    fesa=H1(mesh, order=porder, complex=True, definedon=mesh.Materials("air"))
    for jp in range(0,nphi):
        phi[jp]=(cphi-appphi/2)+appphi*jp/(nphi-1)
        xhat=(np.cos(phi[jp]),np.sin(phi[jp]))
        Eout = exp(-1J*kappa*(x*xhat[0]+y*xhat[1]))
        func1=CoefficientFunction(-1j*kappa*(CoefficientFunction(xhat)*nv)*Eout * u)
        uinf1=Integrate(tuple(func1),mesh,order=porder+1,definedon=
                            mesh.Boundaries("scatterer"))
        vv=GridFunction(fesa)
        vv.vec[:]=0
        vv.Set(Eout,BND,definedon=mesh.Boundaries("scatterer"))
        fvv=CoefficientFunction(grad(vv)*grad(u)-kappa*kappa*vv*u)
        uinf2=Integrate(tuple(fvv),mesh,order=porder+1,definedon=mesh.Materials("air"))
        uinf[jp]=exp(1J*np.pi/4)/np.sqrt(8*np.pi*kappa)*(uinf1+uinf2)
    return(uinf,phi)

def helmsol_penetrable(mesh,porder,ncoef,kappa,incp,farp):
    # Input
    # mesh - the NGSpy mesh
    # porder - order of the FE space
    # ncoef - a coefficient function that gives the value of *n* at points in the mesh (in the PML n=1)
    # kappa - wavenunmber
    # incp - parameters for the incident field
    # farp - parameters for the far field
    #
    # Returns
    # uinf - far field matrix
    # phi - measurement angles
    # theta - incident angles
    #
    fes = H1(mesh, order=porder, complex=True)
    u = fes.TrialFunction()
    v = fes.TestFunction()
    a = BilinearForm(fes)
    a += SymbolicBFI(grad(u)*grad(v) - kappa**2*ncoef*u*v)
    a += SymbolicBFI(-1j*kappa*u*v,definedon=mesh.Boundaries("outerbnd"))
    print('Number of DoFs: ',fes.ndof)
    gfu = GridFunction(fes)
    Draw(gfu,mesh,'us')
    with TaskManager():
        a.Assemble()
        Ainv=a.mat.Inverse()   
    uinf=np.zeros((farp["n"],incp["n"]),dtype=complex)
    theta=np.zeros(incp["n"]);
    center=incp["cent"]
    app=incp["app"]
    for ip in range(0,incp["n"]):
        if ip%10==0:
            print("Done ip = ", ip,' of ',incp["n"])
        if incp["n"]==1:
            theta[0]=0.0
        else:
            theta[ip]=(center-app/2)+app*ip/(incp["n"]-1)
        d=[np.cos(theta[ip]),np.sin(theta[ip])]
        with TaskManager():
            b = LinearForm(fes)
            ui=exp(1J*kappa*(d[0]*x+d[1]*y))
            b += SymbolicLFI(kappa*kappa*(ncoef-1)*ui * v)
            b.Assemble()
            gfu.vec.data =  Ainv * b.vec
            Redraw()
        uinf[:,ip],phi=farf2d_penetrable(gfu,mesh,farp,kappa,porder)
    return(uinf,theta,phi)

def helmsol_dirichlet(mesh,porder,omega,incp,farp):
    fes = H1(mesh, order=porder, complex=True, dirichlet="dirichlet")
    print('Using ',fes.ndof,' DoFs')
    u = fes.TrialFunction()
    v = fes.TestFunction()
    a = BilinearForm(fes)
    a += grad(u)*grad(v)*dx - omega**2*u*v*dx
    a += -1j*omega*u*v*ds("outerbnd")
    
    gfu = GridFunction(fes)
    with TaskManager():
        a.Assemble()
        Ainv=a.mat.Inverse(fes.FreeDofs())
    
    uinf=np.zeros((farp["n"],incp["n"]),dtype=complex)
    center=incp["cent"]
    app=incp["app"]
    if incp["n"]==1:
        theta=np.zeros(1)
        theta[0]=0.0
    else:
        theta=np.linspace(center-app/2,center+app/2,incp["n"])
    for ip in range(0,incp["n"]):
        print('.',end='',flush=True)
        d=[np.cos(theta[ip]),np.sin(theta[ip])]
        ui=exp(1J*omega*(d[0]*x+d[1]*y))
        gfu = GridFunction(fes)
        gfu.vec[:]=0
        gfu.Set(-ui, definedon=mesh.Boundaries("dirichlet"))
        with TaskManager():
            b = LinearForm(fes)
            b.Assemble()
            res = b.vec.CreateVector()
            res.data = - a.mat * gfu.vec
            gfu.vec.data +=  Ainv * res
        if ip%10==0:
            Draw(gfu,mesh,'us')
            Redraw()
        uinf[:,ip],phi=farf2d_dirichlet(gfu,mesh,farp,omega,porder)
    print(' ')
    return(uinf,theta,phi)

def farf2d_dirichlet(u,mesh,farp,omega,porder):
    nphi=farp["n"]
    cphi=farp["cent"]
    appphi=farp["app"]
    if nphi==1:
        phi=np.zeros(1)
        phi[0]=0.0
    else:
        phi=np.linspace(cphi-appphi/2.,cphi+appphi/2,nphi)
    nv=specialcf.normal(mesh.dim)
    uinf=np.zeros(nphi,dtype=complex)
    fesa=H1(mesh, order=porder, complex=True, definedon=mesh.Materials("air"))
    for jp in range(0,nphi):
        xhat=(np.cos(phi[jp]),np.sin(phi[jp]))
        Eout = exp(-1J*omega*(x*xhat[0]+y*xhat[1]))
        func1=CoefficientFunction(-1j*omega*(CoefficientFunction(xhat)*nv)*Eout * u)
        uinf1=Integrate(func1,mesh,order=porder+1,definedon=
                            mesh.Boundaries("dirichlet"))
        vv=GridFunction(fesa)
        vv.vec[:]=0
        vv.Set(Eout,BND,definedon=mesh.Boundaries("dirichlet"))
        fvv=CoefficientFunction(grad(vv)*grad(u)-omega*omega*vv*u)
        uinf2=Integrate(fvv,mesh,order=porder+1,definedon=mesh.Materials("air"))
        uinf[jp]=exp(1J*np.pi/4)/np.sqrt(8*np.pi*omega)*(uinf1+uinf2)
    return(uinf,phi)

