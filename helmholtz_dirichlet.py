# Code to compute multistatic far field pattern for the Dirichlet
# scattering problem for several shapes
from netgen.geom2d import SplineGeometry
from ngsolve import *
from geometries_dirichlet import circle, square, triangle, pac_man, hexagon, Lshape, two_squares
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
ngsglobals.msg_level = 0

from helmsol import farf2d_dirichlet,helmsol_dirichlet



if __name__ == "__main__":
    omega=6  # Wavenumber.  Wavelength of the radiation is 2*pi/omega
    porder=5  # degree of polynmials in the FE method
    NGP=30    # number of triangles per wavelength
    hmax=min((porder*2.*np.pi/omega/NGP,2.*np.pi/omega/5))
    pml_rad=1.5+2*np.pi/omega
    pml_width=2*np.pi/omega/2
    pml_alpha=4.*1j/omega
    print('PML inner radius: ',pml_rad,'  PML width: ',pml_width)
    print('PML alpha: ',pml_alpha)
    print('requested mesh size = ',hmax)
    inc_p={"n":512,"app":2*np.pi,"cent":0} # parameters for incident wave
    far_p={"n":inc_p["n"],"app":inc_p["app"],"cent":0} # params for measurement
    N=ceil(2*omega+1)
    print('Min number of directions for LSM=',N)
    for ngeo in range(0,9):
        if ngeo==0:
            print('Computing for two squares')
            mesh=two_squares(hmax,pml_rad=pml_rad,pml_width=pml_width)
        if ngeo==1:
            print('Computing for Lshape')
            mesh=Lshape(hmax,pml_rad=pml_rad,pml_width=pml_width)
        elif ngeo==2:
            print('Computing for square')
            mesh=square(hmax,pml_rad=pml_rad,pml_width=pml_width)
        elif ngeo==3:
            print('Computing for triangle')
            mesh=triangle(hmax,pml_rad=pml_rad,pml_width=pml_width)
        elif ngeo==4:
            angle=np.pi/8
            print('Computing for pacman, angle=',angle)
            mesh=pac_man(hmax,pml_rad=pml_rad,pml_width=pml_width,order=porder,
                             angle=angle)
        elif ngeo==5:
            angle=2*np.pi/8
            print('Computing for pacman, angle=',angle)
            mesh=pac_man(hmax,pml_rad=pml_rad,pml_width=pml_width,order=porder,
                             angle=angle)
        elif ngeo==6:
            print('Computing for hexagon')
            mesh=hexagon(hmax,pml_rad=pml_rad,pml_width=pml_width)
        elif ngeo==7:
            angle=np.pi/16
            print('Computing for pacman, angle=',angle)
            mesh=pac_man(hmax,pml_rad=pml_rad,pml_width=pml_width,order=porder,
                             angle=angle)
        elif ngeo==8:
            print('Computing for circle')
            mesh=circle(hmax,pml_rad=pml_rad,pml_width=pml_width,order=porder)
             
        mesh.SetPML(pml.Radial(rad=pml_rad,alpha=pml_alpha,origin=(0,0)),
                        "pmlregion")
        uinf,theta,phi=helmsol_dirichlet(mesh,porder,omega,inc_p,far_p)
        if ngeo==0:
            print('saving Two_square')
            np.savez('farff_Two_square.npz',omega,uinf,theta,phi)
            sio.savemat('farff_Two_square.mat',{'omega':omega,'uinf':uinf,
                                          'theta':theta,'phi':phi})
        if ngeo==1:
            print('saving Lshape')
            np.savez('farff_Lshape.npz',omega,uinf,theta,phi)
            sio.savemat('farff_Lshape.mat',{'omega':omega,'uinf':uinf,
                                          'theta':theta,'phi':phi})
        elif ngeo==2:
            print('saving square data')
            np.savez('farf_square.npz',omega,uinf,theta,phi)
            sio.savemat('farf_square.mat',{'omega':omega,'uinf':uinf,
                                               'theta':theta,'phi':phi})
        elif ngeo==3:
            print('saving triangle data')
            np.savez('farf_triangle.npz',omega,uinf,theta,phi)
            sio.savemat('farf_triangle.mat',{'omega':omega,'uinf':uinf,
                                                 'theta':theta,'phi':phi})
        elif ngeo==4:
            print('saving pacman1')
            np.savez('farf_pacman1.npz',omega,uinf,theta,phi,angle)
            sio.savemat('farf_pacman1.mat',{'omega':omega,'uinf':uinf,
                                    'theta':theta,'phi':phi,'angle':angle})
        elif ngeo==5:
            print('saving pacman2')
            np.savez('farf_pacman2.npz',omega,uinf,theta,phi,angle)
            sio.savemat('farf_pacman2.mat',{'omega':omega,'uinf':uinf,
                                    'theta':theta,'phi':phi,'angle':angle})
        elif ngeo==6:
            print('saving hexagon')
            np.savez('farf_hexagon.npz',omega,uinf,theta,phi,angle)
            sio.savemat('farf_hexagon.mat',{'omega':omega,'uinf':uinf,
                                    'theta':theta,'phi':phi})
        elif ngeo==7:
            print('saving pacman3')
            np.savez('farf_pacman3.npz',omega,uinf,theta,phi,angle)
            sio.savemat('farf_pacman3.mat',{'omega':omega,'uinf':uinf,
                                    'theta':theta,'phi':phi,'angle':angle})
        elif ngeo==8:
            print('saving circle data')
            np.savez('farf_circle.npz',omega,uinf,theta,phi)
            sio.savemat('farf_circle.mat',{'omega':omega,'uinf':uinf,
                                               'theta':theta,'phi':phi})

        #plt.figure()
        #plt.contourf(abs(uinf))
        #plt.show()

 
