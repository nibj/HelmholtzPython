# run with "python3 helmholtz_penetrable_test.py"
import sys
import netgen.gui
from netgen.geom2d import SplineGeometry
from ngsolve import *
from geometries_penetrable import circle
import numpy as np
import scipy.special as scs

from helmsol import *

def solve_helmholtz(kappa,ncoef,phi,R,N,npoints,x0):
    #
    # Use special function series
    import cmath as cm
    a=np.zeros(2*N+1,dtype=complex)
    b=np.zeros(2*N+1,dtype=complex)
    rhs=np.zeros(2,dtype=complex)
    for n in range(-N,N+1):
        A=np.zeros((2,2),dtype=complex)
        A[0,0]=scs.hankel1(n,kappa*R)
        A[0,1]=-scs.jv(n,kappa*np.sqrt(ncoef)*R)
        rhs[0]=-(1J**n)*scs.jv(n,kappa*R)*np.exp(-1j*n*phi)
        A[1,0]=kappa*scs.h1vp(n,kappa*R,1)
        A[1,1]=-kappa*np.sqrt(ncoef)*scs.jvp(n,kappa*np.sqrt(ncoef)*R,1)
        rhs[1]=-(1J**n)*kappa*scs.jvp(n,kappa*R)*np.exp(-1j*n*phi)
        ab=np.linalg.solve(A,rhs)
        a[N+n]=ab[0]
        b[N+n]=ab[1]
    theta=np.linspace(0,2*np.pi,npoints)
    d=[np.cos(phi),np.sin(phi)]
    factor=(cm.exp(-1J*np.pi/4)*np.sqrt(2/(np.pi*kappa)))
    factor=factor*np.exp(-1J*kappa*(x0[0]*(np.cos(theta)-d[0])
                                        +x0[1]*(np.sin(theta)-d[1])))
    farff=np.zeros(npoints,dtype=complex)
    for n in range(-N,N+1):
        farff=farff+a[n+N]*np.exp(1J*n*(theta-np.pi/2))
    farff=farff*factor
    return(farff,theta,a,b)

def born_approx(k, R, xlim, ncoef, phi, Ngrid, npoints,x0, ge0):
    # undebugged Born approximation
    # Seems OK in magnitude not phase
    vert_step = 2*xlim/Ngrid # Vertical discretization step size
    hor_step = 2*xlim/Ngrid # Horizontal discretization step size
    Cfac = vert_step*hor_step*np.exp(1J*np.pi/4)*np.sqrt(k**3/(np.pi*8)) # Compute constant factor

    #y1 = np.linspace(x0[0]-xlim, x0[0]+xlim, Ngrid) # Discretize horizontal grid region
    #y2 = np.linspace(x0[1]-xlim, x0[1]+xlim, Ngrid) # Discretize vertical grid region
    y1 = np.linspace(-xlim, xlim, Ngrid) # Discretize horizontal grid region
    y2 = np.linspace(-xlim, xlim, Ngrid) # Discretize vertical grid region
    y1_grid, y2_grid = np.meshgrid(y1, y2) # Produce matrix of coordinate pairs
    y1_flat = y1_grid.flatten() # Flatten into vector
    y2_flat = y2_grid.flatten() # Flatten into vector
    m = np.zeros(Ngrid * Ngrid, dtype=complex) # Initialize refractive index function

    # Solve for refractive index at each grid point
    for i in range(Ngrid * Ngrid):
        if (y1_flat[i]-x0[0])**2 + (y2_flat[i]-x0[1])**2 < R**2:
            m[i] = ncoef - 1
    m_matrix = m.reshape(Ngrid, Ngrid)
    #plt.figure()
    #plt.matshow(np.real(m_matrix))
    #plt.title('True m')
    #plt.show()
    # Angle values to evaluate far field at
    angles = np.linspace(0, 2 * np.pi, npoints)
    dir_cos = np.cos(angles) # xhat(1)
    dir_sin = np.sin(angles) # xhat(2)
    d_cos_phi = np.cos(phi) # d(1)
    d_sin_phi = np.sin(phi) # d(2)

    # Exponential matrix
    A = np.zeros((npoints, Ngrid * Ngrid), dtype=complex)
    for i in range(npoints):
        xhat_dot_y = (dir_cos[i] - d_cos_phi) * y1_flat + (dir_sin[i] - d_sin_phi) * y2_flat
        A[i, :] = np.exp(1j * k * xhat_dot_y)
    # Far field pattern
    A = Cfac * A
    farff = A @ m

    return farff, angles, m, A

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    kappa=8# Wavenumber
    porder=4 # order of the polynmials in the FEM
    nval=1.1 # n inside scatterer.  Born works for small nval
    x0=(0.0,0.0)  # FEM and Born should work for off-centered circles
    R=1
    
    inc_p={"n":1,"app":2*np.pi,"cent":0} # parameters for incident wave
    far_p={"n":129,"app":inc_p["app"],"cent":0} # params for measurement
    
    hmax_s=2*np.pi/kappa/np.sqrt(nval)/8 # mesh size in scatterer
    hmax_a=2*np.pi/kappa/8 # mesh size in air
    #
    print('hmax_s=',hmax_s)
    print('hmax_a=',hmax_a)
    print('wavelength lambda=',2*np.pi/kappa)
    geo="circle"
    
    pml_rad=1+2*2*np.pi/kappa  # radius of inner edge of the PML
    pml_delta=2*np.pi/kappa  # thickness of the PML
    pml_parameter=1j
    mesh=circle(hmax_s=hmax_s,hmax_a=hmax_a,pml_rad=pml_rad,
                            pml_delta=pml_delta,order=porder,x0=x0,R=R)

    print('Materials present: ',mesh.GetMaterials())
    print('Boundary labels: ',mesh.GetBoundaries())
    mesh.SetPML(pml.Radial(rad=pml_rad,alpha=pml_parameter,origin=(0,0)),
                    "pmlregion")
    ncoef=CoefficientFunction([nval if mat=='D' else 1
                              for mat in mesh.GetMaterials()])
    scatter=CoefficientFunction([0 if mat=="D" else 1 if
                                mat=="air" else 2 
                                for mat in mesh.GetMaterials()])
    #Draw(scatter,mesh,'scatterer') # use to check the gtometry is OK
   
    uinf,theta,phi=helmsol_penetrable(mesh,porder,ncoef,kappa,inc_p,far_p) # FEM
    phi=0
    N=100
    npoints=far_p["n"]
    uinfs,theta,a,b=solve_helmholtz(kappa,nval,phi,R,N,npoints,x0) # series
    xlim=4
    Ngrid=200
    uinf_born, theta_born , m , A = born_approx(kappa, R, xlim, nval, phi, Ngrid, npoints,x0, geo) # Born
    plt.figure()
    plt.plot(np.real(uinf),label='FEM real part')
    plt.plot(np.imag(uinf),label='FEM imag part')
    plt.plot(np.real(uinfs),label='Series real part')
    plt.plot(np.imag(uinfs),label='Series imag part')
    plt.plot(np.real(uinf_born),label='Born real part')
    plt.plot(np.imag(uinf_born),label='Born imag part')
    plt.legend()
    plt.show()

    plt.figure()
    plt.plot(np.abs(uinf),label='abs(FEM)')
    plt.plot(np.abs(uinfs),label='abs(Series)')
    plt.plot(np.abs(uinf_born),label='abs(Born)')
    plt.legend()
    plt.show()

    plt.figure()
    plt.plot(np.abs(np.squeeze(uinf[:,0]-uinfs)),label='err(FEM-series)')
    plt.legend()
    plt.title('Error abs(FE-Series)')
    plt.show()
  
