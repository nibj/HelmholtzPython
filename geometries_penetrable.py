from netgen.geom2d import SplineGeometry
from ngsolve import *
import netgen.gui
import numpy as np


def circle(hmax_s=0.3,hmax_a=0.3,pml_rad=2,R=1,pml_delta=.8,order=3,x0=(0,0)):
    geo = SplineGeometry()
    geo.AddCircle( x0, pml_rad+pml_delta, leftdomain=2, bc="outerbnd")
    geo.AddCircle( x0, pml_rad, leftdomain=1, rightdomain=2, bc="innerbnd")
    geo.AddCircle( x0,R,leftdomain=3,rightdomain=1, bc="scatterer")
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    geo.SetMaterial(3,'D')
    geo.SetDomainMaxH(3,hmax_s)
    geo.SetDomainMaxH(2,hmax_a)
    geo.SetDomainMaxH(1,hmax_a)
    mesh = Mesh(geo.GenerateMesh ())
    mesh.Curve(order)
    return(mesh)

def peanut(hmax_s=0.6,hmax_a=0.6,pml_rad=5,pml_delta=.8,order=3):
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_delta, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    Curve= lambda t: (1.2*sqrt(3*cos(t*2*np.pi)**2+1)*cos(t*2*np.pi),
                        1.2*sqrt(3*cos(t*2*np.pi)**2+1)*sin(t*2*np.pi))
    geo.AddCurve(Curve,leftdomain=3,rightdomain=1,bc="scatterer")
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    geo.SetMaterial(3,'D')
    geo.SetDomainMaxH(3,hmax_s)
    geo.SetDomainMaxH(2,hmax_a)
    geo.SetDomainMaxH(1,hmax_a)
    mesh = Mesh(geo.GenerateMesh ())
    mesh.Curve(order)
    return(mesh)

def square(hmax_s=0.1,hmax_a=0.1,pml_rad=2,pml_delta=.8,order=3):
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_delta, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    
    # start code for square
    L2=1/2
    p1,p2,p3,p4 = [ geo.AppendPoint(x,y) for x,y in [ (-L2,-L2), (L2,-L2),
                                                           (L2,L2),(-L2,L2)]]
    geo.Append (["line", p1, p2],leftdomain=3,rightdomain=1,bc="scatterer")
    geo.Append (["line", p2, p3],leftdomain=3,rightdomain=1,bc="scatterer")
    geo.Append (["line", p3, p4],leftdomain=3,rightdomain=1,bc="scatterer")
    geo.Append (["line", p4, p1],leftdomain=3,rightdomain=1,bc="scatterer")
    geo.SetMaterial(3,'D')
    geo.SetDomainMaxH(3,hmax_s)
    # end code for square
    
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    geo.SetDomainMaxH(1,hmax_a)
    geo.SetDomainMaxH(2,hmax_a)
    # genenerate an NGSpy mesh
    mesh = Mesh(geo.GenerateMesh(maxh=hmax_a))
    mesh.Curve(order) # The PML boundaries are curved.
    return(mesh)

def ellipses_in_circle(R_circ=1.5,hmax_a=0.2,hmax_s=.1,pml_delta=.5,
                           pml_rad=1.8,order=3,
            ellip_data={'numellip':3,
        'R1a':0.5,'R1b':.2,'xcen1':(0.3,0.2),'ang1':np.pi/8,'ind1':3,
                             'R2a':0.5,'R2b':.2,'xcen2':(-0.5,-0.2),'ang2':np.pi/8,'ind2':4,
                             'R3a':0.5,'R3b':.2,'xcen3':(0.4,-0.4),'ang3':-np.pi/8,'ind3':5}):
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_delta, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    geo.AddCircle( (0,0), R_circ, leftdomain=6, rightdomain=1,bc="scatterer")
    # first ellipse
    ang=ellip_data["ang1"]
    R1=ellip_data["R1a"]
    R2=ellip_data["R1b"]
    xcen=ellip_data["xcen1"]
    ind=ellip_data["ind1"]
    geo=add_ellipse(geo,R1,R2,xcen,ang,ind,6)
    geo.SetMaterial(ind,"ellip1")
    geo.SetDomainMaxH(ind,hmax_s)
    if ellip_data["numellip"] > 1:
    # second ellipse
        ang=ellip_data["ang2"]
        R1=ellip_data["R2a"]
        R2=ellip_data["R2b"]
        xcen=ellip_data["xcen2"]
        ind=ellip_data["ind2"]
        geo=add_ellipse(geo,R1,R2,xcen,ang,ind,6)
        geo.SetMaterial(ind,"ellip2")
        geo.SetDomainMaxH(ind,hmax_s)
    if ellip_data["numellip"] > 2:
    # third ellipse
        ang=ellip_data["ang3"]
        R1=ellip_data["R3a"]
        R2=ellip_data["R3b"]
        xcen=ellip_data["xcen3"]
        ind=ellip_data["ind3"]
        geo=add_ellipse(geo,R1,R2,xcen,ang,ind,6)
        geo.SetMaterial(ind,"ellip3")
        geo.SetDomainMaxH(5,hmax_s)
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    geo.SetMaterial(6,"big_circle")
    geo.SetDomainMaxH(6,hmax_s)
    geo.SetDomainMaxH(2,hmax_a)
    geo.SetDomainMaxH(1,hmax_a)
    mesh = Mesh(geo.GenerateMesh())
    mesh.Curve(order)
    return(mesh)


                             
def ellipses(hmax_a=0.1,hmax_s=.8,pml_delta=.5,pml_rad=1.8,order=3,
                 ellip_data={'numellip':3,
                             'R1a':0.5,'R1b':.2,'xcen1':(0.3,0.2),'ang1':np.pi/8,'ind1':3,
                             'R2a':0.5,'R2b':.2,'xcen2':(-0.5,-0.2),'ang2':np.pi/8,'ind2':4,
                             'R3a':0.5,'R3b':.2,'xcen3':(0.4,-0.4),'ang3':-np.pi/8,'ind3':5}):
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_delta, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    # first ellipse
    ang=ellip_data["ang1"]
    R1=ellip_data["R1a"]
    R2=ellip_data["R1b"]
    xcen=ellip_data["xcen1"]
    ind=ellip_data["ind1"]
    geo=add_ellipse(geo,R1,R2,xcen,ang,ind,1)
    geo.SetMaterial(ind,"ellip1")
    geo.SetDomainMaxH(ind,hmax_s)
    if ellip_data["numellip"] > 1:
    # second ellipse
        ang=ellip_data["ang2"]
        R1=ellip_data["R2a"]
        R2=ellip_data["R2b"]
        xcen=ellip_data["xcen2"]
        ind=ellip_data["ind2"]
        geo=add_ellipse(geo,R1,R2,xcen,ang,ind,1)
        geo.SetMaterial(ind,"ellip2")
        geo.SetDomainMaxH(ind,hmax_s)
    if ellip_data["numellip"] > 2:
    # third ellipse
        ang=ellip_data["ang3"]
        R1=ellip_data["R3a"]
        R2=ellip_data["R3b"]
        xcen=ellip_data["xcen3"]
        ind=ellip_data["ind3"]
        geo=add_ellipse(geo,R1,R2,xcen,ang,ind,1)
        geo.SetMaterial(ind,"ellip3")
        geo.SetDomainMaxH(5,hmax_s)
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    geo.SetDomainMaxH(2,hmax_a)
    geo.SetDomainMaxH(1,hmax_a)
    mesh = Mesh(geo.GenerateMesh())
    mesh.Curve(order)
    return(mesh)

def add_ellipse(geo,R1,R2,xcen,ang,ind,ind_out):
    # Here the curve defines an ellipse.  note t isin [0,1]
    mat=np.array([[np.cos(ang),np.sin(ang)],[-np.sin(ang),np.cos(ang)]])
    print('ind',ind,'ang=',ang,'x0',xcen,'R1',R1,'R2',R2)
    Curve= lambda t: (
        xcen[0]+mat[0,0]*R1*cos(t*2*np.pi)+mat[0,1]*R2*sin(t*2*np.pi),
        xcen[1]+mat[1,0]*R1*cos(t*2*np.pi)+mat[1,1]*R2*sin(t*2*np.pi))
    if ind_out==1:
        geo.AddCurve(Curve,leftdomain=ind,rightdomain=ind_out,bc="scatterer")
    else:
        geo.AddCurve(Curve,leftdomain=ind,rightdomain=ind_out)
    return(geo)

def random_scatterer(hmax=0.1,pml_rad=2.,pml_width=.8,xcen=(0,0)):
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_width, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    Rmax=1
    Rmin=.2
    npoints=np.random.randint(3,high=15)
    print('Number of points= ',npoints)
    angs=np.random.rand(npoints)
    angs=2*np.pi*np.sort(angs)
    mags=np.random.rand(npoints)
    Rads=(Rmax-Rmin)*mags+Rmin
    xav=0
    yav=0
    for j in range(0,npoints):
        xav=xav+Rads[j]*np.cos(angs[j])
        yav=yav+Rads[j]*np.sin(angs[j])
    # used to place origin inside the scatterer
    xav=xav/npoints
    yav=yav/npoints
    for j in range(0,npoints):
        x=Rads[j]*np.cos(angs[j])-xav
        y=Rads[j]*np.sin(angs[j])-yav
        Rads[j]=np.sqrt(x*x+y*y)
        angs[j]=np.arctan2(x,y)
    Ra=np.zeros((npoints,2))
    Ra[:,0]=angs
    Ra[:,1]=Rads
    Ras=np.sort(Ra,0)
    Rads=Ras[:,1]
    angs=Ras[:,0]
    print('Rads = ',Rads)
    print('angs = ',angs)
    for j in range(0,npoints-1):
        x1=Rads[j]*np.cos(angs[j])+xcen[0]
        y1=Rads[j]*np.sin(angs[j])+xcen[1]
        x2=Rads[j+1]*np.cos(angs[j+1])+xcen[0]
        y2=Rads[j+1]*np.sin(angs[j+1])+xcen[1]
        p1=geo.AppendPoint(x1,y1)
        p2=geo.AppendPoint(x2,y2)
        geo.Append (["line", p1, p2],leftdomain=3,rightdomain=1,bc="scatterer",
                    maxh=hmax/2)
    x1=Rads[npoints-1]*np.cos(angs[npoints-1])+xcen[0]
    y1=Rads[npoints-1]*np.sin(angs[npoints-1])+xcen[1]
    x2=Rads[0]*np.cos(angs[0])+xcen[0]
    y2=Rads[0]*np.sin(angs[0])+xcen[1]
    p1=geo.AppendPoint(x1,y1)
    p2=geo.AppendPoint(x2,y2)
    geo.Append (["line", p1, p2],leftdomain=3,rightdomain=1,bc="scatterer",
                    maxh=hmax/2)
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    geo.SetMaterial(3, "D")
    mesh = Mesh(geo.GenerateMesh(maxh=hmax))
    return(mesh,angs,Rads)

def Lshape(hmax_a=0.1,pml_rad=2.,pml_width=.8,L2=np.sqrt(np.pi)/2,xcen=(0,0),
               ang=0):
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_width, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    #L2=np.sqrt(np.pi)/2 # for equal area
    mat=np.array([[np.cos(ang),np.sin(ang)],[-np.sin(ang),np.cos(ang)]])
    X1=np.matmul(mat,np.array([-L2,-L2]))+xcen
    X2=np.matmul(mat,np.array([-L2,L2]))+xcen
    X3=np.matmul(mat,np.array([L2,L2]))+xcen
    X4=np.matmul(mat,np.array([L2,0]))+xcen
    X5=np.matmul(mat,np.array([0,0]))+xcen
    X6=np.matmul(mat,np.array([0,-L2]))+xcen
    p1,p2,p3,p4,p5,p6 = [ geo.AppendPoint(x,y) for x,y in
                    [ X1, X2,X3, X4,X5,X6]]
    geo.Append (["line", p2, p1],leftdomain=3,rightdomain=1,bc="scatterer",
                    maxh=hmax/2)
    geo.Append (["line", p3, p2],leftdomain=3,rightdomain=1,bc="scatterer",
                    maxh=hmax/2)
    geo.Append (["line", p4, p3],leftdomain=3,rightdomain=1,bc="scatterer",
                    maxh=hmax/2)
    geo.Append (["line", p5, p4],leftdomain=3,rightdomain=1,bc="scatterer",
                    maxh=hmax/2)
    geo.Append (["line", p6, p5],leftdomain=3,rightdomain=1,bc="scatterer",
                    maxh=hmax/2)
    geo.Append (["line", p1, p6],leftdomain=3,rightdomain=1,bc="scatterer",
                    maxh=hmax/2)
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    geo.SetMaterial(3, "D")
    geo.SetDomainMaxH(3,hmax_s)
    geo.SetDomainMaxH(2,hmax_a)
    geo.SetDomainMaxH(1,hmax_a)
    mesh = Mesh(geo.GenerateMesh(maxh=hmax_a))
    return(mesh)

def Ushape(hmax_a=0.1,hmax_s=0.05,pml_rad=2.,pml_width=.8,
               Box=[-.5,.3,-.3,.3],t=0.1, xcen=(0,0),ang=0):
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_width, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    #
    mat=np.array([[np.cos(ang),np.sin(ang)],[-np.sin(ang),np.cos(ang)]])
    Lmin=Box[0]
    Lmax=Box[1]
    Hmin=Box[2]
    Hmax=Box[3]
    print(Lmin,Lmax,Hmin,Hmax)
    X1=np.matmul(mat,np.array([Lmin,Hmin]))+xcen
    X2=np.matmul(mat,np.array([Lmax,Hmin]))+xcen
    X3=np.matmul(mat,np.array([Lmax,Hmin+t]))+xcen
    X4=np.matmul(mat,np.array([Lmin+t,Hmin+t]))+xcen
    X5=np.matmul(mat,np.array([Lmin+t,Hmax-t]))+xcen
    X6=np.matmul(mat,np.array([Lmax,Hmax-t]))+xcen
    X7=np.matmul(mat,np.array([Lmax,Hmax]))+xcen
    X8=np.matmul(mat,np.array([Lmin,Hmax]))+xcen
    p1,p2,p3,p4,p5,p6,p7,p8 = [ geo.AppendPoint(x,y) for x,y in
                    [ X1, X2,X3, X4,X5,X6,X7,X8]]
    geo.Append (["line", p1, p2],leftdomain=3,rightdomain=1,bc="scatterer")
    geo.Append (["line", p2, p3],leftdomain=3,rightdomain=1,bc="scatterer")
    geo.Append (["line", p3, p4],leftdomain=3,rightdomain=1,bc="scatterer")
    geo.Append (["line", p4, p5],leftdomain=3,rightdomain=1,bc="scatterer")
    geo.Append (["line", p5, p6],leftdomain=3,rightdomain=1,bc="scatterer")
    geo.Append (["line", p6, p7],leftdomain=3,rightdomain=1,bc="scatterer")
    geo.Append (["line", p7, p8],leftdomain=3,rightdomain=1,bc="scatterer")
    geo.Append (["line", p8, p1],leftdomain=3,rightdomain=1,bc="scatterer")
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    geo.SetMaterial(3, "D")
    geo.SetDomainMaxH(3,hmax_s)
    geo.SetDomainMaxH(2,hmax_a)
    geo.SetDomainMaxH(1,hmax_a)
    mesh = Mesh(geo.GenerateMesh(maxh=hmax_a))
    return(mesh)

def Venn(hmax=0.1,pml_rad=3.,pml_width=.8,x1=0.8,R1=1.1,order=3):
# intersecting circles centered at (-x1,R1) and (x1,R1) and radius R1
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_width, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    #
    p1=geo.AppendPoint(0,np.sqrt(R1**2-x1**2))
    print(x1+R1,-x1+R1)
    alpha=(x1+R1)/np.sqrt(R1**2-x1**2)
    p1a=geo.AppendPoint(x1+R1,np.sqrt(R1**2-x1**2)+alpha*x1)
    p1aR=geo.AppendPoint(-(x1+R1),np.sqrt(R1**2-x1**2)+alpha*x1)
    p2=geo.AppendPoint(x1+R1,0)
    p2R=geo.AppendPoint(-(x1+R1),0)
    p2a=geo.AppendPoint(x1+R1,-(np.sqrt(R1**2-x1**2)+alpha*x1))
    p2aR=geo.AppendPoint(-(x1+R1),-(np.sqrt(R1**2-x1**2)+alpha*x1))
    p3=geo.AppendPoint(0,-np.sqrt(R1**2-x1**2))
    beta=(-x1+R1)/np.sqrt(R1**2-x1**2)
    p3a=geo.AppendPoint(-x1+R1,np.sqrt(R1**2-x1**2)-beta*x1)
    p3aR=geo.AppendPoint(x1-R1,np.sqrt(R1**2-x1**2)-beta*x1)
    print(-x1+R1,np.sqrt(R1**2-x1**2)-beta*x1)
    p4=geo.AppendPoint(-x1+R1,0)
    p4a=geo.AppendPoint(-x1+R1,-(np.sqrt(R1**2-x1**2)-beta*x1))
    p4aR=geo.AppendPoint(x1-R1,-(np.sqrt(R1**2-x1**2)-beta*x1))
    p4R=geo.AppendPoint(x1-R1,0)
    geo.Append(["spline3",p2,p1a,p1],leftdomain=3,rightdomain=1,bc="scatterer")
    geo.Append(["spline3",p3,p2a,p2],leftdomain=3,rightdomain=1,bc="scatterer")
    geo.Append(["spline3",p1,p3a,p4],leftdomain=3,rightdomain=4)
    geo.Append(["spline3",p4,p4a,p3],leftdomain=3,rightdomain=4)
    geo.Append(["spline3",p3,p4aR,p4R],leftdomain=5,rightdomain=4)
    geo.Append(["spline3",p4R,p3aR,p1],leftdomain=5,rightdomain=4)
    geo.Append(["spline3",p1,p1aR,p2R],leftdomain=5,rightdomain=1,bc="scatterer")
    geo.Append(["spline3",p2R,p2aR,p3],leftdomain=5,rightdomain=1,bc="scatterer")
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    geo.SetMaterial(3, "DL")
    geo.SetMaterial(4, "DLR")
    geo.SetMaterial(5, "DR")
    mesh=Mesh(geo.GenerateMesh(maxh=hmax))
    mesh.Curve(order)
    return(mesh)
    
   
def Shepp_Logan(hmax=0.05,pml_rad=1.,pml_width=.8,order=3):
    # Modified Shepp-Logan type phantom (no subdomains overlap)
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_width, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    # a
    xcen=(0,0)
    ind=3
    ind_out=1
    ang=0
    R1=.69
    R2=.92
    add_ellipse(geo,R1,R2,xcen,ang,ind,ind_out)
    #b
    xcen=(0,-0.0184)
    ind=4
    ind_out=3
    ang=0
    R1=.6624
    R2=.874
    add_ellipse(geo,R1,R2,xcen,ang,ind,ind_out)
    #c
    xcen=(0.22,0)
    ind=5
    ind_out=4
    ang=18./360.*2.*np.pi
    R1=.11
    R2=.31
    add_ellipse(geo,R1,R2,xcen,ang,ind,ind_out)
    #d
    xcen=(-0.22,0)
    ind=6
    ind_out=4
    ang=-18./360.*2.*np.pi
    R1=.16
    R2=.41
    add_ellipse(geo,R1,R2,xcen,ang,ind,ind_out)
    #e
    xcen=(0,0.45)
    ind=7
    ind_out=4
    ang=0
    R1=.21
    R2=.25
    add_ellipse(geo,R1,R2,xcen,ang,ind,ind_out)
    #f
    xcen=(0,0.1)
    ind=8
    ind_out=4
    ang=0
    R1=.046
    R2=.046
    add_ellipse(geo,R1,R2,xcen,ang,ind,ind_out)
    #g
    xcen=(0,-0.1+.1)
    ind=9
    ind_out=4
    ang=0
    R1=.046
    R2=.046
    add_ellipse(geo,R1,R2,xcen,ang,ind,ind_out)
    #h
    xcen=(-0.08,-0.605)
    ind=10
    ind_out=4
    ang=0
    R1=.046
    R2=.023
    add_ellipse(geo,R1,R2,xcen,ang,ind,ind_out)
    #i
    xcen=(0.0,-0.605)
    ind=11
    ind_out=4
    ang=0
    R1=.023
    R2=.023
    add_ellipse(geo,R1,R2,xcen,ang,ind,ind_out)
    #j
    xcen=(0.06,-0.605)
    ind=12
    ind_out=4
    ang=0
    R1=.023
    R2=.046
    add_ellipse(geo,R1,R2,xcen,ang,ind,ind_out)
    #
    #
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    geo.SetMaterial(3, "a")
    geo.SetMaterial(4, "b")
    geo.SetMaterial(5, "c")
    geo.SetMaterial(6, "d")
    geo.SetMaterial(7, "e")
    geo.SetMaterial(8, "f")
    geo.SetMaterial(9, "g")
    geo.SetMaterial(10, "h")
    geo.SetMaterial(11, "i")
    geo.SetMaterial(12, "j")
    mesh=Mesh(geo.GenerateMesh(maxh=hmax))
    mesh.Curve(order)
    return(mesh)

def Bullseye(hmax=0.1,pml_rad=2.,pml_width=.8,order=3):
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_width, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    # a
    xcen=(0,0)
    ind=3
    ind_out=1
    ang=0
    R1=1.
    R2=1.
    add_ellipse(geo,R1,R2,xcen,ang,ind,ind_out)
    #b
    xcen=(0,0)
    ind=4
    ind_out=3
    ang=0
    R1=.8
    R2=.8
    add_ellipse(geo,R1,R2,xcen,ang,ind,ind_out)
    #c
    xcen=(0,-0.0184)
    ind=5
    ind_out=4
    ang=0
    R1=.6
    R2=.4
    add_ellipse(geo,R1,R2,xcen,ang,ind,ind_out)
    #
    #
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    geo.SetMaterial(3, "a")
    geo.SetMaterial(4, "b")
    geo.SetMaterial(5, "c")
    mesh=Mesh(geo.GenerateMesh(maxh=hmax))
    mesh.Curve(order)
    return(mesh)

def Squares(hmax=0.1,pml_rad=2.,pml_width=.8):
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_width, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    # points
    # need x1<x2<x3
    x1=-.4
    x2=.4
    x3=.6
    # need y1<y2<y3<y4
    y1=-.4
    y2=-.2
    y3=.1
    y4=.4
    p1=geo.AppendPoint(0,y2)
    p2=geo.AppendPoint(0,y3)
    p3=geo.AppendPoint(x2,y3)
    p4=geo.AppendPoint(x2,y2)
    # a
    geo.Append(["line",p1,p2],leftdomain=3,rightdomain=4)
    geo.Append(["line",p2,p3],leftdomain=3,rightdomain=4)
    geo.Append(["line",p3,p4],leftdomain=5,rightdomain=4)
    geo.Append(["line",p4,p1],leftdomain=3,rightdomain=4)
    #
    # b
    p5=geo.AppendPoint(x3,y3)
    p6=geo.AppendPoint(x3,y2)
    geo.Append(["line",p5,p3],leftdomain=5,rightdomain=1,bc="scatterer")
    geo.Append(["line",p6,p5],leftdomain=5,rightdomain=1,bc="scatterer")
    geo.Append(["line",p4,p6],leftdomain=5,rightdomain=1,bc="scatterer")
    #c
    p7=geo.AppendPoint(x2,y1)
    p8=geo.AppendPoint(x2,y4)
    p9=geo.AppendPoint(x1,y4)
    p10=geo.AppendPoint(x1,y1)
    geo.Append(["line",p3,p8],leftdomain=3,rightdomain=1,bc="scatterer")
    geo.Append(["line",p8,p9],leftdomain=3,rightdomain=1,bc="scatterer")
    geo.Append(["line",p9,p10],leftdomain=3,rightdomain=1,bc="scatterer")
    geo.Append(["line",p10,p7],leftdomain=3,rightdomain=1,bc="scatterer")
    geo.Append(["line",p7,p4],leftdomain=3,rightdomain=1,bc="scatterer")
    #
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    geo.SetMaterial(3, "a")
    geo.SetMaterial(4, "b")
    geo.SetMaterial(5, "c")
    mesh=Mesh(geo.GenerateMesh(maxh=hmax))
    return(mesh)

# Use this code to test each mesh type
if __name__ == "__main__":
    print('meshing Venn')
    #mesh=ellipses_in_circle()
    #mesh,angs,Rads=random_scatterer()
    #mesh=Lshape()
    #mesh=Ushape()
    #mesh=ellipses()
    #geom="Shepp"
    #geom="Venn"
    #geom="Bullseye"
    geom="Squares"
    if geom=="Shepp":
        mesh=Shepp_Logan()
        print(mesh.GetMaterials())
        scatter=CoefficientFunction([2 if mat=="pmlregion" else
                                     1 if mat=="air" else
                                     3 if mat=="a" else
                                     4 if mat=="b" else
                                     6 if mat=="c" else
                                     6 if mat=="d" else
                                     7 if mat=="e" else
                                     8 if mat=="f" else
                                     9 if mat=="g" else
                                     10 if mat=="h" else
                                     11 if mat=="i" else
                                     12 for mat in mesh.GetMaterials()])
    #                                4 if mat=="ellip2" else
    #                                5 if mat=="ellip3" else
    #                                6 if mat=="big_circle" else
    #                                7 for mat in mesh.GetMaterials()])
        V=H1(mesh,order=1)
        gfu=GridFunction(V)
        gfu.Set(1,BND,definedon=mesh.Boundaries("scatterer"))
        Draw(gfu,mesh,'sctt')
    elif geom=="Venn":
        mesh=Venn()
        scatter=CoefficientFunction([2 if mat=="pmlregion" else
                                     1 if mat=="air" else
                                     3 if mat=="DL" else
                                     4 if mat=="DLR" else
                                     5  for mat in mesh.GetMaterials()])
        V=H1(mesh,order=1)
        gfu=GridFunction(V)
        gfu.Set(1,BND,definedon=mesh.Boundaries("scatterer"))
        Draw(gfu,mesh,'sctt')
        nv=specialcf.normal(mesh.dim)
        Draw(nv[1],mesh,'normal')
    elif geom=="Bullseye":
        mesh=Bullseye()
        V=H1(mesh,order=1)
        gfu=GridFunction(V)
        gfu.Set(1,BND,definedon=mesh.Boundaries("scatterer"))
        scatter=CoefficientFunction([2 if mat=="pmlregion" else
                                     1 if mat=="air" else
                                     3 if mat=="a" else
                                     4 if mat=="b" else
                                     5  for mat in mesh.GetMaterials()])
        Draw(gfu,mesh,'sctt')
    else:
        mesh=Squares()
        V=H1(mesh,order=1)
        gfu=GridFunction(V)
        gfu.Set(1,BND,definedon=mesh.Boundaries("scatterer"))
        scatter=CoefficientFunction([2 if mat=="pmlregion" else
                                     1 if mat=="air" else
                                     3 if mat=="a" else
                                     4 if mat=="b" else
                                     5  for mat in mesh.GetMaterials()])

        Draw(gfu,mesh,'sctt')
    Draw(scatter,mesh,'domains')
    input()

