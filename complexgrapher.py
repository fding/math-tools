from PIL import Image
from math import *
import os

screen=(1024,1280)
pix_range=(screen[0]/2,screen[1]/2)
        
# Tabulates values of function on points of a grid
def tabulate(func,rrange,irange,size=(screen[0]/2,screen[1]/2)):
    pix_range=size
    out=[]
    m1=(rrange[1]-rrange[0])/float(pix_range[1])
    m2=(irange[1]-irange[0])/float(pix_range[0])
    for i in range(pix_range[0]):
        row=[]
        for j in range(pix_range[1]):
            real=rrange[0]+m1*j
            imag=irange[0]+m2*i
            try:
                val=func(complex(real,imag))
                if val==float('inf'):
                    row.append([-1,0])
                    continue
            except:
                row.append([-1,0]) #-1 magnitude signals that the function attains the maximum value there, whatever that maybe.
                #We avoid just throwing a large number, because that would ruin the graph.
                continue
            mag=abs(val)
            arg=atan2(val.imag,val.real)
            arg=arg if arg>0 else 2*pi+arg
            row.append([mag,arg])
            
        out.append(row)
    return out
    
# returns color corresponding to modulus, so that color changes continuously
def arg2color(arg):
    if arg<2*pi/3.0:
        n=255*arg/(2*pi/3.0)
        curcolor=[255.0-n,n,0]
    elif arg<4*pi/3.0:
        n=255*(arg-2*pi/3.0)/(2*pi/3.0)
        curcolor=[0,255.0-n,n]
    else:
        n=255*(arg-4*pi/3.0)/(2*pi/3.0)
        curcolor=[n,0,255.0-n]
    return curcolor

def plot_table(nums,size=(screen[0]/2,screen[1]/2),mode='',brightness=0.75):
    pix_range=size
    if mode=='':
        cgraph=Image.new('RGBA',pix_range,(0,0,0))
    else: 
        cgraph=Image.new(mode,pix_range,(0,0,0))
    max_r=0
    avg_mag=0.0
    for i in range(pix_range[0]):
        for j in range(pix_range[1]):
            mag=nums[i][j][0]
            max_r=max_r if max_r>mag else mag
            avg_mag+=mag
            
    avg_mag=avg_mag/(pix_range[0]*pix_range[1])
    grid=cgraph.load()
    gamma=1 # Scaling factor for intensity of colors.
    try:
        gamma=log(brightness)/log(avg_mag/max_r)
    except:
        pass
    print gamma
    for i in range(pix_range[0]):
        for j in range(pix_range[1]):
            mag=nums[i][j][0]
            if mag<0:
                mag=max_r #Makes the singularities attain maximum value.
            arg=nums[i][j][1]
            m=(mag*1.0/max_r)**gamma
            color=arg2color(arg)
            color=[int(color[0]),int(color[1]),int(color[2])]
            if mode=='RGB':
                color=[int(color[0]*m),int(color[1]*m),int(color[2]*m)]
            if mode=='RGBA':
                color.append(int(255*(1-m)))
            if mode=='':
                color=[int(color[0]*m),int(color[1]*m),int(color[2]*m)]
                color.append(int(255*(1-m)))
            grid[i,j]=tuple(color)
    cgraph=cgraph.rotate(90)
    cgraph.save('temp_from_complexgrapher.png')
    return cgraph

def plot_function(func,rrange,irange,size=(screen[0]/2,screen[1]/2),mode='',brightness=0.75):
    '''
        Plots a complex function f. Each complex number is given as a magnitude and an argument.
        The output is a PIL Image of the graph. Each point z is colored
        according to the magnitude and argument of f(z): the color indicates the argument, with pure red being
        real (pure green and pure blue are along the other cube roots of unity). The magnitude is shown by the intensity:
        Very bright regions have large magnitudes, and dark regions having small magnitudes.

        Brightness is an argument that gives the brightness (out of 1) of the mediumly bright points.
    '''
    return plot_table(tabulate(func,rrange,irange,size),size,mode,brightness)
    
def square_plot(func,rrange,irange,rorange,iorange,hlines=7,vlines=7,size=(screen[1]/2,screen[0]/2),resolution=1):
    '''
        Plots a square with hlines horizontal lines and vlines vertical lines in the input plane,
        and outputs the image of that square under func.
        Resolution, an integer, determines the density of input points.
    '''
    cgraph=Image.new('1',size,True)
    domain=[]
    mr=(rrange[1]-rrange[0])/float(size[0]*resolution)
    mi=(irange[1]-irange[0])/float(hlines)
    for i in range(hlines):
        imag=mi*i+irange[0]
        for j in range(resolution*size[0]):
            real=mr*j+rrange[0]
            domain.append(complex(real,imag))
    
    mi=(irange[1]-irange[0])/float(size[1]*resolution)
    mr=(rrange[1]-rrange[0])/float(vlines)
    for i in range(vlines):
        real=mr*i+rrange[0]
        for j in range(resolution*size[1]):
            imag=mi*j+irange[0]
            domain.append(complex(real,imag))
            
    output=[]
    for i in domain:
        try:
            output.append(func(i))
        except:
            output.append(complex(rorange[1]+1,iorange[1]+1))
    out=cgraph.load()
    mr=(size[0]-1)/(rorange[1]-rorange[0])
    mi=(size[1]-1)/(iorange[1]-iorange[0])
    for i in output:
        if i.real>rorange[1] or i.real<rorange[0]:
            continue
        if i.imag>iorange[1] or i.imag<iorange[0]:

            continue
        x=(i.real-rorange[0])*mr
        y=(i.imag-iorange[0])*mi
        out[int(x),int(y)]=False
    
    cgraph.save('temp_from_complexgrapher.png')
    os.startfile('temp_from_complexgrapher.png')
    inputgraph=Image.new('1',size,True)
    points=inputgraph.load()
    for i in domain:
        points[int((size[0]-1)/(rrange[1]-rrange[0])*(i.real-rrange[0])),int((size[1]-1)/(irange[1]-irange[0])*(i.imag-irange[0]))]=False
    
    return inputgraph,cgraph
