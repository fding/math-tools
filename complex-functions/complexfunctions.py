'''
    A module for defining various complex functions
'''
from math import *
import cmath

i=complex(0,1)

def Theta(z,tau):
    '''
        The standard Jacoby Theta function. Tau should have an imaginary part
        greater than zero.
    '''
    if tau.imag<=0:
        raise
    q=cmath.exp(pi*i*tau)
    zeta=cmath.exp(2*pi*i*z)
    cur=complex(1,0)
    cur1=complex(1,0)
    threshold=1e-7
    maxrecursion=10000
    n=1
    curq=complex(1,0)
    curzeta1=complex(1,0)
    curzeta2=curzeta1
    while True:
        cur1=cur
        curq=q**(n*n)
        curzeta1*=zeta
        curzeta2/=zeta
        n+=1
        cur+=curq*(curzeta1+curzeta2)
        if abs(cur-cur1)<threshold*abs(cur):
            return cur
        if n>maxrecursion:
            return float('inf')
        
        
def Theta01(z,tau):
    return Theta(z+0.5,tau)
def Theta10(z,tau):
    return cmath.exp(0.25*pi*i*tau+pi*i*z)*Theta(z+0.5*tau,tau)
def Theta11(z,tau):
    return cmath.exp(0.25*pi*i*tau+pi*i*(z+0.5))*Theta(z+0.5*tau+0.5,tau)
    
def w(z,tau):
    #A handy tool to calculate WElliptic
    if tau.imag<=0:
        raise
    #Theta functions rise quickly with large imaginary parts.
    #We know that w is periodic, so we make z smaller.
    while -z.imag>tau.imag:
        z+=tau
    while z.imag>tau.imag:
        z-=tau
    pi2=pi*pi
    t1=Theta(0,tau)
    t2=Theta10(0,tau)
    t1*=t1
    t2*=t2
    t3=Theta01(z,tau)
    t4=Theta11(z,tau)
    t5=t3/t4
    t5*=t5
    return pi2*t1*t2*t5-pi2/3.0*(t1*t1+t2*t2)
    
def Zeta(s):
    '''
        The Riemann zeta function
    '''
    gamma=[.57721566490153286,-0.072815845483676724860586,
           -0.0096903631928723184845303,0.002053834420303345866160,
            0.0023253700654673000574,0.0007933238173010627017,
           -0.00023876934543019960986,-0.0005272895670577510,
           -0.00035212335380,-0.0000343947744,0.000205332814909,
           .27018443954e-3,.16727291210514019,-.20920926205929994e-4,
           -.283468655320241e-3,-.1996968583089697747e-3,.26277037109918e-4,
           .307368408149e-3,.503605453047355629e-3,.4663435615115594e-3] #20 values
    
    if s==0: return -0.5
    if s==1: return float('inf')
    if s.real<0:
        return 2**s*pi**(s-1)*cmath.sin(pi*s/2)*Gamma(1-s)*Zeta(1-s)
    if s.real>0 and s.real<1:
        if s.real<0.5:
            return 2**s*pi**(s-1)*cmath.sin(pi*s/2)*Gamma(1-s)*Zeta(1-s)
        cursum=0.0
        n=1
        curdem=1
        curnum=1
        while True:
            curnum*=(s+n-1)
            curdem*=n+1
            curterm=(Zeta(s+n)-1)*curnum/float(curdem)
            cursum+=curterm
            n+=1
            if abs(curterm)<1e-4:
                return s/(s-1.0)-cursum
    cursum=0
    n=1
    while True:
        curTerm=(-1)**(n-1)/float(n)**s
        cursum+=curTerm
        n+=1
        if abs(curTerm)<1e-8:
            return 1.0/(1.0-2**(1-s))*cursum

def Stieltjes(n):
    m=100
    out=0
    for k in range(1,m+1):
        out+=log(k)**n/float(k)
    return out-log(m)**(n+1)/float(n+1)

def Gamma(z):
    '''
        Gamma function.
    '''
    if z.imag==0:
        if z.real<0:
            if int(z.real)==z.real:
                return float('inf')
    if z.real>1:
        return (z-1)*Gamma(z-1)
    if z.real<0:
        return Gamma(z+1.0)/z
    p=[1.000000000190015,76.18009172947146,-86.50532032941677,24.01409824083091,
       -1.231739572450155,1.208650973866179e-3,-5.395239384953e-6]
    cursum=p[0]
    for n in range(1,7):
        cursum+=p[n]/(z+float(n))
    return (sqrt(2*pi)/z)*cursum*(z+5.5)**(z+0.5)*cmath.exp(-(z+5.5))


def WElliptic(z,w1,w2):
    '''
        The Weierstrauss Elliptic function with periods w1 and w2. w1/w2 should be non-real.
    '''
    if (w2/w1).imag<=0:
        w1,w2=w2,w1
    return w(z/w1,w2/w1)/w1/w1


'''
    The next few functions are the Jacobi elliptic functions
'''
def tauf(k):
    '''
        An auxilliary function to find tau from k
    '''
    k=cmath.sqrt(1-k*k)
    k=cmath.sqrt(k)
    l=0.5*(1.0-k)/(1.0+k)
    q=l+2.0*l**5+15.0*l**9+150.0*l**13+1707.*l**17+20910.*l**21+268616.*l**25
    return cmath.log(q)/(i*pi)

    
def sn(u,k):
    '''
        k is the elliptic modulus
    '''
    tau=tauf(k)
    t=Theta(0,tau)
    z=u/pi/t**2
    return -t*Theta11(z,tau)/Theta10(0,tau)/Theta01(z,tau)
    
def cn(u,k):
    tau=tauf(k)
    t=Theta(0,tau)
    z=u/pi/t**2
    return Theta01(0,tau)*Theta10(z,tau)/Theta10(0,tau)/Theta01(z,tau)
    
def dn(u,k):
    tau=tauf(k)
    t=Theta(0,tau)
    z=u/pi/t**2
    return Theta01(0,tau)*Theta(z,tau)/Theta(0,tau)/Theta01(z,tau)
    
