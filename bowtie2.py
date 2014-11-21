#!/usr/bin/python

#Open Source Bow Tie
#Create a pdf containing a pattern for Bow Tie
#Published under GNU GPL v3
#Crappy python/postscript coding by
#Juho Happola, Mari Ijas
#juho.happola@iki.fi

#requires: numpy, ps2pdf to convert the newly created ps into pdf

#caveats: this has not been tested with invalid input

import datetime
import numpy
import pylab
import os
import sys
import scipy
import copy
import math
import random

band = 2.5
height = 8.0
width = 12.0
joinseam = 3.0
knot = 4.0
necksize = 38.0
seam = 0.7
fabricthickness = 0.3
wrinkle = 0.5
fcoefs = [1, 0.15, 0.02, 0, 0, 0]
version = '0.2'

tilt = random.gauss(-0.75,0.1)*math.pi/2
tilt2 = random.gauss(-1.25,0.1)*math.pi/2


seam = seam+fabricthickness

def readSpecs(filename):
  global width
  global height
  global necksize
  global band
  global knot
  global joinseam
  global seam
  global fcoefs
  global wrinkle
  f = open(filename,'r')
  lines = f.readlines()
  efl = []
  for line in lines:
    if '#' in line:
      efl.append(line[0:line.find('#')])
    else:
      efl.append(line[0:-1])
  width = float(efl[0])
  height = float(efl[1])
  necksize = float(efl[2])
  band = float(efl[3])
  knot = float(efl[4])
  joinseam = float(efl[5])
  seam = float(efl[6])
  fabricthickness = float(efl[7])
  wrinkle = float(efl[8])
  seam = seam + fabricthickness
  fcoefs = [float(efl[9])]
  if (len(efl))>10:
    for num in range(10,len(efl)):
      fcoefs.append(float(efl[num]))
  f.close()

if (len(sys.argv) > 1):
  readSpecs(sys.argv[1])

N = int(1e2)
tilt = random.gauss(-0.75,0.1)*math.pi/2
tilt2 = random.gauss(-1.25,0.1)*math.pi/2

A = 0.25*(height-band)
D = width/0.6666667
L = 2.0
C = 0.5*band
AC = 0.5*(height-band-wrinkle)

fcoefs = pylab.array(fcoefs)*(1.0/(sum(fcoefs)))

B = 0.5*necksize-1.0*L+0.5*knot+joinseam

k=3*math.pi/D

def psHeader():
  return '%%!PS-Adobe-2.0\n%%%%BoundingBox: 0 0 595 842\n%%%%DocumentMedia: 210x294mm 595 842 0 () ()\n\n/m { moveto } bind def\n/l { lineto } bind def\n/s { stroke } bind def\n1 setlinewidth\n0.1 setgray\n\ngsave\n'

def centerPolygon(poly,amp):
  xave = 0.5*(max(poly[0])+min(poly[0]))
  yave = 0.5*(max(poly[1])+min(poly[1]))
  xl = (max(poly[0])-min(poly[0]))
  yl = (max(poly[1])-min(poly[1]))
  xpap = 21.0
  ypap = 29.4
  xtol = xpap-xl
  ytol = ypap-yl
  return (-1.0*xave+0.5*xpap+random.gauss(0,amp*xtol), -1.0*yave+0.5*ypap+random.gauss(0,amp*ytol))

def psFooter():
  return 'grestore \n /Times findfont 14 scalefont setfont 400 50 moveto ( Open Source Bow Tie %s %d) show \n showpage\n%%%%EOF'%(version,datetime.datetime.now().year)

def psLine(x1,y1,x2,y2):
  xx1=x1*(595.0/21.0)
  xx2=x2*(595.0/21.0)
  yy1=y1*(848.0/29.4)
  yy2=y2*(848.0/29.4)
  return '%f %f m %f %f l s\n\n'%(xx1,yy1,xx2,yy2)

def drawPolygon(f,poly,shift):
  ss = poly.shape
  ss = ss[1]
  for index in range(0,ss-1):
    xx1 = poly[0][index]+shift[0]
    yy1 = poly[1][index]+shift[1]
    xx2 = poly[0][index+1]+shift[0]
    yy2 = poly[1][index+1]+shift[1]
    f.write(psLine(xx1,yy1,xx2,yy2))
  xx1 = poly[0][ss-1]+shift[0]
  yy1 = poly[1][ss-1]+shift[1]
  xx2 = poly[0][0]+shift[0]
  yy2 = poly[1][0]+shift[1]
  f.write(psLine(xx1,yy1,xx2,yy2))

def printBowtie(shapes):
  filu = open('bowtie.ps','w')
  filu.write(psHeader())
  corr=centerPolygon(shapes[0],0.1)
  drawPolygon(filu,shapes[0],corr)
  drawPolygon(filu,shapes[1],corr)
  corr=centerPolygon(shapes[2],0.1)
  drawPolygon(filu,shapes[2],corr)
  drawPolygon(filu,shapes[3],corr)  
  filu.write(psFooter())
  filu.close()

def polygonHull(poly,dist):
  ss = poly.shape
  ss = ss[1]

  xs = []
  ys = []
  xx0 = poly[0][0]
  yy0 = poly[1][0]
  xx1 = poly[0][ss-1]
  yy1 = poly[1][ss-1]
  xx2 = poly[0][1]
  yy2 = poly[1][1]
  dy = xx2-xx1
  dx = yy1-yy2
  direction = pylab.array([[dx],[dy]])
  cc = dist/pylab.norm(direction)
  xs.append(xx0+cc*direction[0])
  ys.append(yy0+cc*direction[1])

  for index in range(1,ss-1):
    xx0 = poly[0][index]
    yy0 = poly[1][index]
    xx1 = poly[0][index-1]
    yy1 = poly[1][index-1]
    xx2 = poly[0][index+1]
    yy2 = poly[1][index+1]
    dy = xx2-xx1
    dx = yy1-yy2
    direction = pylab.array([[dx],[dy]])
    cc = dist/pylab.norm(direction)
    xs.append(xx0+cc*direction[0])
    ys.append(yy0+cc*direction[1])

  xx0 = poly[0][ss-1]
  yy0 = poly[1][ss-1]
  xx1 = poly[0][ss-2]
  yy1 = poly[1][ss-2]
  xx2 = poly[0][0]
  yy2 = poly[1][0]
  dy = xx2-xx1
  dx = yy1-yy2
  direction = pylab.array([[dx],[dy]])
  cc = dist/pylab.norm(direction)
  xs.append(xx0+cc*direction[0])
  ys.append(yy0+cc*direction[1])

  xv = pylab.array(xs)
  yv = pylab.array(ys)
  rv1 = pylab.reshape(xv,(1,len(xv)))
  rv2 = pylab.reshape(yv,(1,len(yv)))
  return pylab.concatenate((copy.deepcopy(rv1),copy.deepcopy(rv2)),0)

rot = pylab.array([[math.cos(tilt), -1.0*math.sin(tilt)], [math.sin(tilt), math.cos(tilt)]])
rot2 = pylab.array([[math.cos(tilt2), -1.0*math.sin(tilt2)], [math.sin(tilt2), math.cos(tilt2)]])
dang = pylab.array([[0.0, -1.0], [1.0, 0]])
tm = pylab.zeros([N,N])

for ind in range(0,N):
  tm[N-ind-1][ind]=1.0

x1=pylab.linspace(-L,0,N)
x2=pylab.linspace(0,(1.0/3.0)*D,N)
x8=pylab.linspace((1.0/3.0)*D,D,N)
x3=D*pylab.ones(N)
x4=-1.0*L*pylab.ones(N)

x5=pylab.linspace(0,B,N)
x6=B*pylab.ones(N)
x7=0.0*pylab.ones(N)

#yh = pylab.concatenate((A*pylab.ones(N/3),AC*pylab.ones(2*N/3)))
y1 = C*pylab.ones(N)
y2 = pylab.zeros(N)
for l in range(0,len(fcoefs)):
  y2 = y2-1.0*(fcoefs[l]*pylab.cos((l+1)*k*x2))
y2 = y2-y2[0]
y2 = (2.0*A/max(y2))*y2
y2 = y2+C
y7 = pylab.zeros(N)
for l in range(0,len(fcoefs)):
  y7 = y7-1.0*(fcoefs[l]*pylab.cos((l+1)*k*x8))
y7 = y7-min(y7)
y7 = (1.0/max(y7))*y7
y7 = (y2[-1]-C-wrinkle/2.0)*y7
y7 = y7 + C+wrinkle/2.0
y3 = pylab.linspace(y7[-1],-1.0*y7[-1],N)
y4 = -1.0*y2;
y5 = -1.0*y1;
y6 = pylab.linspace(-C,C,N);

xx = pylab.concatenate((copy.deepcopy(x1),copy.deepcopy(x2),copy.deepcopy(x8),copy.deepcopy(x3),copy.deepcopy(pylab.dot(tm,x8)),copy.deepcopy(pylab.dot(tm,x2)), copy.deepcopy(pylab.dot(tm,x1)), copy.deepcopy(x4)),0)
yy = pylab.concatenate((copy.deepcopy(y1),copy.deepcopy(y2),copy.deepcopy(y7),copy.deepcopy(y3),copy.deepcopy(-1.0*pylab.dot(tm,y7)),copy.deepcopy(pylab.dot(tm,y4)),copy.deepcopy(pylab.dot(tm,y5)),copy.deepcopy(y6)),0)

xneck = pylab.concatenate((copy.deepcopy(x5),copy.deepcopy(x6),copy.deepcopy(pylab.dot(tm,x5)),copy.deepcopy(x7)),0)
yneck = pylab.concatenate((copy.deepcopy(y1),copy.deepcopy(pylab.dot(tm,y6)),copy.deepcopy(y5),copy.deepcopy(y6)),0)

xxn = pylab.reshape(xx,(1,len(xx)))
yyn = pylab.reshape(yy,(1,len(yy)))

xxneck = pylab.reshape(xneck,(1,len(xneck)))
yyneck = pylab.reshape(yneck,(1,len(yneck)))

tht = pylab.concatenate((xxn,yyn),0)
points = pylab.dot(rot,tht)
shapes = []
shapes.append(copy.deepcopy(points))

tht = polygonHull(points,seam)
shapes.append(copy.deepcopy(tht))

tht = pylab.concatenate((xxneck,yyneck),0)
points=pylab.dot(rot2,tht)
shapes.append(copy.deepcopy(points))

tht = polygonHull(points,seam)
shapes.append(copy.deepcopy(tht))

printBowtie(shapes)
os.system('ps2pdf bowtie.ps')