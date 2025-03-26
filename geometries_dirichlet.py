import netgen.gui
from netgen.geom2d import SplineGeometry
from ngsolve import *
import numpy as np
## need ellipse and allow rotation for all shapes.
def circle(hmax=0.1,pml_rad=2.,pml_width=.8,order=3,R=1,xcen=(0,0)):
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_width, leftdomain=2, rightdomain=0,
                       bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    geo.AddCircle( xcen ,R,leftdomain=0,rightdomain=1, bc="dirichlet",
                       maxh=hmax/2)
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    mesh = Mesh(geo.GenerateMesh (maxh=hmax))
    mesh.Curve(order)
    return(mesh)

def hexagon(hmax=0.1,pml_rad=2.,pml_width=.8,L2=1,xcen=(0,0),ang=0):
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_width, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    #L2=1 # for equal area
    #L2=1
    c=cos(np.pi/3.)
    s=sin(np.pi/3)
    mat=np.array([[np.cos(ang),np.sin(ang)],[-np.sin(ang),np.cos(ang)]])
    X1=np.matmul(mat,np.array([-L2,0]))+xcen
    X2=np.matmul(mat,np.array([-L2*c,L2*s]))+xcen
    X3=np.matmul(mat,np.array([L2*c,L2*s]))+xcen
    X4=np.matmul(mat,np.array([L2,0]))+xcen
    X5=np.matmul(mat,np.array([L2*c,-L2*s]))+xcen
    X6=np.matmul(mat,np.array([-L2*c,-L2*s]))+xcen
    p4a=(L2+xcen[0],xcen[1])
    p5a=(L2*c+xcen[0],-L2*s+xcen[1])
    p6a=(-L2*c+xcen[0],-L2*s+xcen[1])
    p1,p2,p3,p4,p5,p6 = [ geo.AppendPoint(x,y) for x,y in [ X1,X2,X3,X4,X5,X6]]
    geo.Append (["line", p2, p1],leftdomain=0,rightdomain=1,bc="dirichlet",
                    maxh=hmax/2)
    geo.Append (["line", p3, p2],leftdomain=0,rightdomain=1,bc="dirichlet",
                    maxh=hmax/2)
    geo.Append (["line", p4, p3],leftdomain=0,rightdomain=1,bc="dirichlet",
                    maxh=hmax/2)
    geo.Append (["line", p5, p4],leftdomain=0,rightdomain=1,bc="dirichlet",
                    maxh=hmax/2)
    geo.Append (["line", p6, p5],leftdomain=0,rightdomain=1,bc="dirichlet",
                    maxh=hmax/2)
    geo.Append (["line", p1, p6],leftdomain=0,rightdomain=1,bc="dirichlet",
                    maxh=hmax/2)
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    mesh = Mesh(geo.GenerateMesh(maxh=hmax))
    return(mesh)

def random_scatterer(hmax=0.1,pml_rad=2.,pml_width=.8,xcen=(0,0)):
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_width, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    Rmax=1
    Rmin=.2# for equal area
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
        geo.Append (["line", p1, p2],leftdomain=0,rightdomain=1,bc="dirichlet",
                    maxh=hmax/2)
    x1=Rads[npoints-1]*np.cos(angs[npoints-1])+xcen[0]
    y1=Rads[npoints-1]*np.sin(angs[npoints-1])+xcen[1]
    x2=Rads[0]*np.cos(angs[0])+xcen[0]
    y2=Rads[0]*np.sin(angs[0])+xcen[1]
    p1=geo.AppendPoint(x1,y1)
    p2=geo.AppendPoint(x2,y2)
    geo.Append (["line", p1, p2],leftdomain=0,rightdomain=1,bc="dirichlet",
                    maxh=hmax/2)
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    mesh = Mesh(geo.GenerateMesh(maxh=hmax))
    return(mesh,angs,Rads)


def ellipse(hmax=0.1,pml_rad=2.,pml_width=.8,R1=1,R2=1,xcen=(0,0),ang=0):
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_width, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    #
    mat=np.array([[np.cos(ang),np.sin(ang)],[-np.sin(ang),np.cos(ang)]])
    X1=np.matmul(mat,np.array([-R1,0]))+xcen
    X2=np.matmul(mat,np.array([-R1,R2]))+xcen
    X3=np.matmul(mat,np.array([0,R2]))+xcen
    X4=np.matmul(mat,np.array([R1,R2]))+xcen
    X5=np.matmul(mat,np.array([R1,0]))+xcen
    X6=np.matmul(mat,np.array([R1,-R2]))+xcen
    X7=np.matmul(mat,np.array([0,-R2]))+xcen
    X8=np.matmul(mat,np.array([-R1,-R2]))+xcen
    p1,p2,p3,p4,p5,p6,p7,p8 = [ geo.AppendPoint(x,y) for x,y in [ X1,X2,
                                            X3,X4,X5,X6,X7,X8]]
    geo.Append (["spline3", p1, p2, p3],leftdomain=1,rightdomain=0,bc='dirichlet')
    geo.Append (["spline3", p3, p4, p5],leftdomain=1,rightdomain=0,bc='dirichlet')
    geo.Append (["spline3", p5, p6, p7],leftdomain=1,rightdomain=0,bc='dirichlet')
    geo.Append (["spline3", p7, p8, p1],leftdomain=1,rightdomain=0,bc='dirichlet')
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    mesh = Mesh(geo.GenerateMesh(maxh=hmax))
    return(mesh)

def square(hmax=0.1,pml_rad=2.,pml_width=.8,L2=np.sqrt(np.pi)/2,xcen=(0,0),ang=0):
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_width, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    #L2=np.sqrt(np.pi)/2 # for equal area
    #L2=1
    mat=np.array([[np.cos(ang),np.sin(ang)],[-np.sin(ang),np.cos(ang)]])
    X1=np.matmul(mat,np.array([-L2,-L2]))+xcen
    X2=np.matmul(mat,np.array([L2,-L2]))+xcen
    X3=np.matmul(mat,np.array([L2,L2]))+xcen
    X4=np.matmul(mat,np.array([-L2,L2]))+xcen
    p1,p2,p3,p4 = [ geo.AppendPoint(x,y) for x,y in [ X1,X2,X3,X4]]
    geo.Append (["line", p1, p2],leftdomain=0,rightdomain=1,bc="dirichlet",
                    maxh=hmax/2)
    geo.Append (["line", p2, p3],leftdomain=0,rightdomain=1,bc="dirichlet",
                    maxh=hmax/2)
    geo.Append (["line", p3, p4],leftdomain=0,rightdomain=1,bc="dirichlet",
                    maxh=hmax/2)
    geo.Append (["line", p4, p1],leftdomain=0,rightdomain=1,bc="dirichlet",
                    maxh=hmax/2)
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    mesh = Mesh(geo.GenerateMesh(maxh=hmax))
    return(mesh)

def Lshape(hmax=0.1,pml_rad=2.,pml_width=.8,L2=np.sqrt(np.pi)/2,xcen=(0,0),
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
    geo.Append (["line", p2, p1],leftdomain=0,rightdomain=1,bc="dirichlet",
                    maxh=hmax/2)
    geo.Append (["line", p3, p2],leftdomain=0,rightdomain=1,bc="dirichlet",
                    maxh=hmax/2)
    geo.Append (["line", p4, p3],leftdomain=0,rightdomain=1,bc="dirichlet",
                    maxh=hmax/2)
    geo.Append (["line", p5, p4],leftdomain=0,rightdomain=1,bc="dirichlet",
                    maxh=hmax/2)
    geo.Append (["line", p6, p5],leftdomain=0,rightdomain=1,bc="dirichlet",
                    maxh=hmax/2)
    geo.Append (["line", p1, p6],leftdomain=0,rightdomain=1,bc="dirichlet",
                    maxh=hmax/2)
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    mesh = Mesh(geo.GenerateMesh(maxh=hmax))
    return(mesh)


def triangle(hmax=0.1,pml_rad=2.,pml_width=.8,R=1,xcen=(0,0),
                 ang=0,angle=np.pi/8):
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_width, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    mat=np.array([[np.cos(ang),np.sin(ang)],[-np.sin(ang),np.cos(ang)]])
    X1=np.matmul(mat,np.array([-R,0]))+xcen;
    X2=np.matmul(mat,np.array([R*np.cos(angle),R*np.sin(angle)]))+xcen;
    X3=np.matmul(mat,np.array([R*np.cos(angle),-R*np.sin(angle)]))+xcen;
    p1,p2,p3 = [ geo.AppendPoint(x,y) for x,y in [X1, X2, X3]]
    geo.Append (["line", p1, p3],leftdomain=0,rightdomain=1,bc="dirichlet",
                    maxh=hmax/2)
    geo.Append (["line", p3, p2],leftdomain=0,rightdomain=1,bc="dirichlet",
                    maxh=hmax/2)
    geo.Append (["line", p2, p1],leftdomain=0,rightdomain=1,bc="dirichlet",
                    maxh=hmax/2)
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    mesh = Mesh(geo.GenerateMesh(maxh=hmax))
    return(mesh)

def pac_man(hmax=0.1,pml_rad=2.,pml_width=.8,order=3,angle=np.pi/8,
                xcen=np.array([0,0]),
                ang=0,R=1):
    # assume angle< pi/2
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_width, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    # now starts the scatterer
    mat=np.array([[np.cos(ang),np.sin(ang)],[-np.sin(ang),np.cos(ang)]])
    x1=R*np.cos(angle/2)
    y1=R*np.sin(angle/2)
    X1=np.matmul(mat,np.array([x1,y1]))+xcen;
    X2=np.matmul(mat,np.array([R*np.cos(angle/2),-R*np.sin(angle/2)]))+xcen;
    X3=np.matmul(mat,np.array([(R**2-y1**2)/x1,R]))+xcen;
    X4=np.matmul(mat,np.array([(R**2-y1**2)/x1,-R]))+xcen;
    print(X1,X2,X3,X4)
    p1,p2,p3,p4 = [ geo.AppendPoint(x,y) for x,y in [ X1, X2, X3, X4]]
    X5=np.matmul(mat,np.array([0,-R]))+xcen;
    X6=np.matmul(mat,np.array([-R,0]))+xcen;
    X7=np.matmul(mat,np.array([0,R]))+xcen;
    p5,p6,p7 = [ geo.AppendPoint(x,y) for x,y in [ X5, X6, X7 ]]
    X8=np.matmul(mat,np.array([-R,R]))+xcen;
    X9=np.matmul(mat,np.array([-R,-R]))+xcen;
    X10=np.matmul(mat,np.array([0,0]))+xcen;
    p8,p9,p10 = [ geo.AppendPoint(x,y) for x,y in [ X8, X9, X10 ]]
    geo.Append (["line", p10, p1],leftdomain=0,rightdomain=1,
                    bc="dirichlet",maxh=hmax/2)
    geo.Append (["line", p2, p10],leftdomain=0,rightdomain=1,
                    bc="dirichlet",maxh=hmax/2)
    geo.Append (["spline3", p1, p3, p7],leftdomain=0,rightdomain=1,
                    bc="dirichlet",maxh=hmax/2)
    geo.Append (["spline3", p7, p8, p6],leftdomain=0,rightdomain=1,
                    bc="dirichlet",maxh=hmax/2)
    geo.Append (["spline3", p6, p9, p5],leftdomain=0,rightdomain=1,
                    bc="dirichlet",maxh=hmax/2)
    geo.Append (["spline3", p5, p4, p2],leftdomain=0,rightdomain=1,
                    bc="dirichlet",maxh=hmax/2)
    #
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    mesh = Mesh(geo.GenerateMesh(maxh=hmax))
    return(mesh)

def two_squares(hmax=0.1,pml_rad=2.,pml_width=.8, L2=np.sqrt(np.pi)/2,
                    xcen=(0,0), ang=np.pi/4):
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_width, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    #L2=np.sqrt(np.pi)/2 # for equal area
    #L2=1
    DL1=L2
    mat=np.array([[np.cos(ang),np.sin(ang)],[-np.sin(ang),np.cos(ang)]])
    p1a=np.matmul(mat,np.array([-L2,-L2]))+xcen;
    p2a=np.matmul(mat,np.array([-L2+DL1,-L2]))+xcen;
    p3a=np.matmul(mat,np.array([-L2+DL1,-L2+DL1]))+xcen;
    p4a=np.matmul(mat,np.array([-L2,-L2+DL1]))+xcen;
    print('p1a',p1a,' p2a',p2a)
    p1,p2,p3,p4 = [ geo.AppendPoint(x,y) for x,y in [ p1a, p2a,  p3a,p4a]]
    geo.Append (["line", p1, p2],leftdomain=0,rightdomain=1,bc="dirichlet",maxh=hmax/2)
    geo.Append (["line", p2, p3],leftdomain=0,rightdomain=1,bc="dirichlet",maxh=hmax/2)
    geo.Append (["line", p3, p4],leftdomain=0,rightdomain=1,bc="dirichlet",maxh=hmax/2)
    geo.Append (["line", p4, p1],leftdomain=0,rightdomain=1,bc="dirichlet",maxh=hmax/2)

    DL2=L2/1.5
    p1a=np.matmul(mat,np.array([L2,L2]))+xcen;
    p2a=np.matmul(mat,np.array([L2-DL2,L2]))+xcen;
    p3a=np.matmul(mat,np.array([L2-DL2,L2-DL2]))+xcen;
    p4a=np.matmul(mat,np.array([L2,L2-DL2]))+xcen;
 
    p1,p2,p3,p4 = [ geo.AppendPoint(x,y) for x,y in [ p1a,p2a,p3a,p4a]]
    geo.Append (["line", p1, p2],leftdomain=0,rightdomain=1,bc="dirichlet",maxh=hmax/2)
    geo.Append (["line", p2, p3],leftdomain=0,rightdomain=1,bc="dirichlet",maxh=hmax/2)
    geo.Append (["line", p3, p4],leftdomain=0,rightdomain=1,bc="dirichlet",maxh=hmax/2)
    geo.Append (["line", p4, p1],leftdomain=0,rightdomain=1,bc="dirichlet",maxh=hmax/2)
     
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    mesh = Mesh(geo.GenerateMesh(maxh=hmax))
    return(mesh)

# Use this code to test each mesh type
if __name__ == "__main__":
    seed=1
    np.random.seed(seed)    
    j=0
    if j==0:
        for j in range(0,10):
            r_points=np.random.uniform(size=5)
            R=r_points[0]+.5
            xcen=np.array([0.25*(2*r_points[1]-1),0.25*(2*r_points[2]-1)])
            ang=(2*r_points[3]-1)*2*np.pi
            #mesh=circle(R=R,xcen=xcen)
            #mesh=hexagon(L2=R,ang=ang)
            #mesh,angs,rads=random_scatterer(xcen=xcen)
            #mesh=square(L2=R*2/3,xcen=xcen,ang=ang)
            R1=R
            R2=r_points[4]+.5
            #mesh=ellipse(R1=R1,R2=R2,ang=ang)
            #mesh=Lshape(L2=np.sqrt(np.pi)/2+0.3*(2*r_points[4]-1),ang=ang)
            print('R=',R,' xcen=',xcen,' ang=',ang)
            #mesh=two_squares(ang=ang)
            #mesh=pac_man(angle=4*np.pi/16,ang=ang,R=R,xcen=xcen)
            mesh=triangle(R=R,ang=ang,xcen=xcen)
            scatter=CoefficientFunction([0 if mat=="pmlregion" else
                                     1 if mat=="air" else 2 
                              for mat in mesh.GetMaterials()])
            Draw(scatter,mesh,'domains')
            print('mesh boundaries: ',mesh.GetBoundaries())
            print('mesh materials: ',mesh.GetMaterials())
            input()

