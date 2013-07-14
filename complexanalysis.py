from cmath import *

# numerical path integration of a function f along a path.
def PathIntegral(f,path,t0,t1,accuracy=0.00001):
    msum=0
    numrecurse=0
    N=4000
    error=0
    curestimate=0
    while True:
        delta=(t1-t0)/float(N)
        msum=0
        for i in range(N):
            curz=path(delta*0.5+i*delta)
            msum+=(path(i*delta+delta)-path(i*delta))*f(curz)
        if numrecurse>0:
            if abs(msum-curestimate)<accuracy*abs(msum) or abs(msum-curestimate)<accuracy:
                return msum
        if numrecurse>20:
            print msum
            print curestimate
            break
        curestimate=msum
        N+=4000
        numrecurse+=1

# uses Cauchy's formula to compute winding numbers
def WindingNumber(path,t0,t1,z0=0):
    integral=PathIntegral(lambda z: 1.0/(z-z0),path,t0,t1)
    ans=integral/(2.0*pi)
    return int(round(ans.imag))

# uses path integrals to compute Laurent series coefficients
def LaurentCoefficient(f,z0,m,R=1.0):
    r=R*0.5
    integral=PathIntegral(lambda z: f(z)*(z-z0)**(-1-m),lambda t: z0+r*e**(1.0j*t),0.0,2*pi)
    ans=integral/(2.0*pi)
    return ans.imag

# residue of function
def Residue(f,z0,R=1.0):
    return LaurentCoefficient(f,z0,-1,R)

def DirectionalDerivative(f,z0,direction,h=1e-7):
    return (f(z0+h*direction)-f(z0-h*direction))/(2.0*h*direction)

# Checks Cauchy's theorem to determine analyticity of function in area enclosed by loop
def isAnalytic(f,z0,h=1e-7, tolerance=1e-5):
    d1= DirectionalDerivative(f,z0,1.0,h)
    d2= DirectionalDerivative(f,z0,1.0j,h)
    if abs(d1-d2)<tolerance or abs(d1-d2)<tolerance*abs(d1):
        return True
    return False

# complex derivative
def derivative(f,h=1e-7):
    return lambda z: DirectionalDerivative(f,z,1.0,h)

def SignedArea(path,t0,t1):
    integral=PathIntegral(lambda z: z.imag,path,t0,t1)
    return -integral.real

