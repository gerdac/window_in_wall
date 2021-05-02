import numpy as np
from matplotlib import pyplot
import os
import random as rand
#import csv
#import scipy.sparse as spa
#import scipy.sparse.linalg as splag
#from mpl_toolkits.mplot3d import Axes3D

#task:
task="single"   #"single" or "loopDOE" or "loopGB"

#controls:
withHole=False 
withSpotWave=True
withGate=True
withCurv=True
withOverAngle=True
AM=False
withDrieness=True
gateType="persian"   #"persian" or "halfCircle"

#wall dimensions & discretization
wallWidth = 3.5
wallHeight = 3.2
elSize =0.03

#wall waviness up & down
upAmp=0.05
upFreq=2.5
upPhase=np.pi/3.0

downAmp=0.04
downFreq=2.0
downPhase=0.0

#local waviness
spotX=1.75
spotY=1.0
spotAmp=0.04
spotFreq=3.5
spotPhase=0.0
spotRad=0.5
spotEnd=2.0

#circular window
cenOpenX = 1.0
cenOpenY = 0.0
openRad=0.5

#general shape gate/window
gateSize=0.7
gateX=1.75

#pattern diffuision properties
kPres=1.0
kShear=0.2

#optimization variables
loopIters = 50

#printing variables
layerThick=0.025
increLen=0.02
depth=0.05
nCount=0
elCount=0
frames=1000
nozzleSpeed=0.1

#initialization of global variables

xSize = int(wallWidth / elSize) + 1
ySize = int(wallHeight / elSize) + 1

x=np.zeros ((xSize,ySize))
y=np.zeros ((xSize,ySize))
z=np.zeros ((xSize,ySize))
xi=np.zeros ((xSize,ySize))
eta=np.zeros ((xSize,ySize))

Kx=np.zeros ((xSize*ySize,xSize*ySize))
Ky=np.zeros ((xSize*ySize,xSize*ySize))
rhsX=np.zeros (xSize*ySize)
rhsY=np.zeros (xSize*ySize)
dispX=np.zeros((xSize,ySize))
dispY=np.zeros((xSize,ySize))

curv=np.zeros ((xSize,ySize))
overAngle=np.zeros ((xSize,ySize))

elMissing=np.zeros((xSize-1)*(ySize-1))

eps=0.0000000001

gatePoints = []

############################################# MAIN WALL FUNCTIONS #########################################

def main():
   if task=="single":
      makeAwall(-1)
   elif task=="loopDOE":
      loopDOE()   

def makeAwall(ite):
   meshTheWall()
   setKmat()
   if withGate:
      openGate()
   #if withHole:
   #    openHoleTrans()
   waveTheWall()
   if withSpotWave:
      spotWave()
   #if withCurv:
   #   calGaussCurv()
   #writeCSV()
   writeVTK(ite)

############################################# LINEAR ALGEBRA FUNCTIONS #########################################

def globOf(i,j):
   return (i*ySize+j)

def applyXdir(i,j,val):
   global rhsX
   rhsX=rhsX-Kx[globOf(i,j),:]*val
   Kx[globOf(i,j),:]=0
   Kx[:,globOf(i,j)]=0
   Kx[globOf(i,j),globOf(i,j)]=1.0
   rhsX[globOf(i,j)]=val

def applyYdir(i,j,val):
   global rhsY
   rhsY=rhsY-Ky[globOf(i,j),:]*val
   Ky[globOf(i,j),:]=0
   Ky[:,globOf(i,j)]=0
   Ky[globOf(i,j),globOf(i,j)]=1.0
   rhsY[globOf(i,j)]=val

def setKmat():
  
   if gateType=="persian":
      persianGate()
   elif gateType=="halfCircle":
      cirGate()
   else:
      print("unknown Gate Type")
         
   global rhsX
   n=xSize*ySize
   for i in range(n):
      Kx[i,i]=kPres*4.0
      if (i)>0 and (i%ySize)!=0:
         Kx[i,(i-1)]=Kx[i,(i-1)]-kPres         
      if (i<n-1) and(i%ySize)!=(ySize-1):
         Kx[i,(i+1)]=Kx[i,(i+1)]-kPres
      if (n-i)>ySize:
         Kx[i,i+ySize]=Kx[i,i+ySize]-kPres
      if i>(ySize-1):
         Kx[i,i-ySize]=Kx[i,i-ySize]-kPres

      if (i%ySize)==0:
         Kx[i,(i+1)]=Kx[i,(i+1)]-kPres
      if (i%ySize)==ySize-1:
         Kx[i,(i-1)]=Kx[i,(i-1)]-kPres

   for i in range(n):
      Ky[i,i]=kPres*4.0
      if (i)>0 and (i%ySize)!=0:
         Ky[i,(i-1)]=Ky[i,(i-1)]-kPres         
      if (i<n-1) and(i%ySize)!=(ySize-1):
         Ky[i,(i+1)]=Ky[i,(i+1)]-kPres
      if (n-i)>ySize:
         Ky[i,i+ySize]=Ky[i,i+ySize]-kPres
      if i>(ySize-1):
         Ky[i,i-ySize]=Ky[i,i-ySize]-kPres

      if i<ySize:
         Ky[i,(i+ySize)]=Ky[i,(i+ySize)]-kPres
      if i/ySize>xSize-1:
         Ky[i,(i-ySize)]=Ky[i,(i-ySize)]-kPres

   for p in (gatePoints):
      applyXdir(p[0],p[1],p[2])
      applyYdir(p[0],p[1],p[3])
   
   solX=np.linalg.solve(Kx,rhsX)
   solY=np.linalg.solve(Ky,rhsY)

   for i in range(xSize):
      for j in range(ySize):
         dispX[i][j]=solX[globOf(i,j)]
         dispY[i][j]=solY[globOf(i,j)]

   # pyplot.imshow(dispY)
   # pyplot.imshow(Kx)
   # pyplot.imshow(Kx)
   # print(rhsX)
   # pyplot.imshow(Kx)
   # pyplot.imshow(np.linalg.inv(Kx))
   # pyplot.colorbar()
   # pyplot.show()
   

############################################# GENERAL MATH FUNCTIONS #########################################
         
def cart2pol(cenX,cenY,x, y):
    rho = np.sqrt((x-cenX)**2 + (y-cenY)**2)
    phi = np.arctan2((y-cenY), (x-cenX))
    return(rho, phi)

def pol2cart(cenX, cenY, rho, phi):
    x = rho * np.cos(phi) + cenX
    y = rho * np.sin(phi) + cenY
    return(x, y)

def sinWave(amp,freq,phase,inX):
   return(amp*np.sin(2.0*np.pi*freq*inX+phase))

def blend(inAmin,inAmax,inPoint,inMin,inMax):
   answer=((inMax-inPoint)*inAmin+(inPoint-inMin)*inAmax)/(inMax-inMin)
   answer=max(answer,min(inAmin,inAmax))
   answer=min(answer,max(inAmin,inAmax))
   return (answer)

def fade(start,end,val):
   if val<start:
       return (1.0)
   else:
      if val>end:
          return (0.0)
      else:
         return((1.0+np.cos((np.pi/(end-start))*val+(start*np.pi/(start-end))))/2.0)
         
def projectOnVec(v_start,v,p):
   pp=v_start + np.dot((p-v_start),v)/np.dot(v,v)*v
   return(pp,np.dot(pp-v_start,v))
   
############################################# HELPING FUNCTIONS #########################################

def zFinder(xy):
   i=max(min(int(round(xy[0]/elSize)),xSize-1),1)
   j=max(min(int(round(xy[1]/elSize)),ySize-1),1)
   return(z[i][j])

def indexFinder(xy):
   i=max(min(int(round(abs(xy[0]/elSize))),xSize-1),1)
   j=max(min(int(round(abs(xy[1]/elSize))),ySize-1),1)
   return((i,j))
   
def isInHole(xy):
   answer=False
   dist=np.sqrt((cenOpenX-xy[0])*(cenOpenX-xy[0])+(cenOpenY-xy[1])*(cenOpenY-xy[1]))
   if dist<openRad:
      answer=True
   return(answer)

############################################# WALL OPERATIONS #########################################

def findMissingEls():
   counter=0
   for i in range(xSize-1):
      for j in range(ySize-1):
         elX=x[i][j]+(elSize/2.0)
         elY=y[i][j]+(elSize/2.0)
         dist=np.sqrt((cenOpenX-elX)*(cenOpenX-elX)+(cenOpenY-elY)*(cenOpenY-elY))
         if dist<openRad+(elSize/np.sqrt(2.0)):
            elMissing[counter]=1
         counter=counter+1

def calGaussCurv():
   for i in range(1,xSize-1):
      for j in range(1,ySize-1):
         if (not(withHole) or (np.sqrt((cenOpenX-x[i][j])*(cenOpenX-x[i][j])+(cenOpenY-y[i][j])*(cenOpenY-y[i][j]))>(openRad+elSize*2.0))):
            v1=[x[i][j-1]-x[i][j],y[i][j-1]-y[i][j],z[i][j-1]-z[i][j]]
            v2=[x[i+1][j]-x[i][j],y[i+1][j]-y[i][j],z[i+1][j]-z[i][j]]
            v3=[x[i][j+1]-x[i][j],y[i][j+1]-y[i][j],z[i][j+1]-z[i][j]]
            v4=[x[i-1][j]-x[i][j],y[i-1][j]-y[i][j],z[i-1][j]-z[i][j]]
            normal1=np.cross(v1,v2)
            normal2=np.cross(v2,v3)
            normal3=np.cross(v3,v4)
            normal4=np.cross(v4,v1)
            area1=np.linalg.norm(normal1)
            area2=np.linalg.norm(normal2)
            area3=np.linalg.norm(normal3)
            area4=np.linalg.norm(normal4)
            area=(area1+area2+area3+area4)/4.0
            normal=(normal1/area1)+(normal2/area2)+(normal3/area3)+(normal4/area4)
            teta1=np.arccos(np.dot(v1, v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)))
            teta2=np.arccos(np.dot(v2, v3)/(np.linalg.norm(v2)*np.linalg.norm(v3)))
            teta3=np.arccos(np.dot(v3, v4)/(np.linalg.norm(v3)*np.linalg.norm(v4)))
            teta4=np.arccos(np.dot(v4, v1)/(np.linalg.norm(v4)*np.linalg.norm(v1)))
            curv[i][j]=((2*np.pi)-(teta1+teta2+teta3+teta4))/area
            downUnit=[0,-1.0,0]
            thisAngle=(np.arccos(np.dot(normal, downUnit)/(np.linalg.norm(normal)*np.linalg.norm(downUnit))))
            overAngle[i][j]=180.0*((thisAngle-(np.pi/2.0)))/np.pi


def meshTheWall():
   for i in range(xSize):
      for j in range(ySize):
         x[i][j]=i*elSize
         y[i][j]=j*elSize
         xi[i][j]=x[i][j]
         eta[i][j]=y[i][j]

def openHole():
   for i in range(xSize):
      for j in range(ySize):
         polar=cart2pol(cenOpenX,cenOpenY,x[i][j],y[i][j])
         newRho=(polar[0]+np.sqrt(polar[0]*polar[0]+4*openRad*openRad))/2.0
         (x[i][j],y[i][j])=pol2cart(cenOpenX,cenOpenY,newRho,polar[1])

def openHoleTrans():
   for i in range(xSize):
      for j in range(ySize):
         polar=cart2pol(cenOpenX,cenOpenY,x[i][j],y[i][j])
         if polar[0]>openRad:
             newRho=polar[0]-(openRad*openRad)/polar[0]
             (xi[i][j],eta[i][j])=pol2cart(cenOpenX,cenOpenY,newRho,polar[1])

def waveTheWall():
   for i in range(xSize):
      for j in range(ySize):
         upZ=sinWave(upAmp,upFreq,upPhase,xi[i][j])
         downZ=sinWave(downAmp,downFreq,downPhase,xi[i][j])
         z[i][j]=blend(downZ,upZ,eta[i][j],0,wallHeight)

def spotWave():
   for i in range(xSize):
      for j in range(ySize):
         spotZ=sinWave(spotAmp,spotFreq,spotPhase,xi[i][j])
         spotDist=np.sqrt((xi[i][j]-spotX)*(xi[i][j]-spotX)+(eta[i][j]-spotY)*(eta[i][j]-spotY))
         z[i][j]=z[i][j]+fade(spotRad,spotEnd,spotDist)*spotZ
         

def cirGate():
   gateR=gateSize*1.0
   for i in range (xSize):
      for j in range (ySize):
         dist=np.sqrt((x[i][j]-gateX)*(x[i][j]-gateX)+(y[i][j])*(y[i][j]))
         if abs(dist-gateR)<elSize:
            gatePoints.append([i,j,-(gateX-x[i][j]),y[i][j]])

def persianGate():
   r=gateSize*1.0
   for i in range (xSize):
      for j in range (ySize):
         xx=gateX-x[i][j]
         yy=y[i][j]
         region1=yy<r*(np.sqrt(2.5-np.sqrt(2.0))+np.sqrt(2)/2.0)
         equation1=abs(xx)-r<elSize
         region2=yy>r*(np.sqrt(2.5-np.sqrt(2.0))+np.sqrt(2)/2.0) and yy<r*(np.sqrt(2.5-np.sqrt(2.0))+np.sqrt(2))
         equation2=abs((xx)**2+(yy-r*(np.sqrt(2.5-np.sqrt(2.0))+np.sqrt(2)/2.0))**2)-r*r<elSize*elSize
         region3=yy>r*(np.sqrt(2.5-np.sqrt(2.0))+np.sqrt(2)) and yy<r*(np.sqrt(2.5-np.sqrt(2.0))+np.sqrt(7/2))
         equation3=abs((xx+xx*r*np.sqrt(2)/(2*abs(xx)))**2+(yy-r*(np.sqrt(2.5-np.sqrt(2.0))))**2)-4*r*r<elSize*elSize
         if (region1 and equation1) or (region2 and equation2) or (region3 and equation3):
            gatePoints.append([i,j,-xx,yy])
            
            

def openGate():
   for i in range(xSize):
      for j in range(ySize):
         xi[i][j]=x[i][j]-dispX[i][j]
         eta[i][j]=y[i][j]-dispY[i][j]
         #targX=x[i][j]-dispX[i][j]
         #targY=y[i][j]-dispY[i][j]
         #xi[i][j]=x[indexFinder((targX,targY))]
         #eta[i][j]=y[indexFinder((targX,targY))]


############################################# PRINTING PATH #########################################

def AMtopol():
   pnt = open("AMpath.pnt", "w")
   elm = open("AMpath.elm", "w")
   #runners:
   b1l_start=np.array([0,0])
   b1l_end=np.array([0,cenOpenY])
   b1r_start=np.array([wallWidth,0])
   b1r_end=np.array([wallWidth,cenOpenY])

   b2l_start=np.array([0,cenOpenY])
   b2l_end=np.array([0,wallHeight])
   b2r_start=np.array([cenOpenX,cenOpenY])
   b2r_end=np.array([cenOpenX,cenOpenY+openRad])

   b3l_start=np.array([cenOpenX,cenOpenY])
   b3l_end=np.array([cenOpenX,cenOpenY+openRad])
   b3r_start=np.array([wallWidth,cenOpenY])
   b3r_end=np.array([wallWidth,wallHeight])

   b4l_start=np.array([cenOpenX,cenOpenY+openRad])
   b4l_end=np.array([0,wallHeight])
   b4r_start=np.array([cenOpenX,cenOpenY+openRad])
   b4r_end=np.array([wallWidth,wallHeight])

   AMcellWriter(b1l_start,b1l_end,b1r_start,b1r_end,pnt,elm)
   AMcellWriter(b2l_start,b2l_end,b2r_start,b2r_end,pnt,elm)
   AMcellWriter(b3l_start,b3l_end,b3r_start,b3r_end,pnt,elm)
   AMcellWriter(b4l_start,b4l_end,b4r_start,b4r_end,pnt,elm)

   elm.close()
   pnt.close()

def AMcellWriter(v1_start,v1_end,v2_start,v2_end,pnt,elm):

   global nCount
   global elCount
   v1=v1_end-v1_start
   v2=v2_end-v2_start
   leftDiv=int(np.linalg.norm(v1)/layerThick)
   rightDiv=int(np.linalg.norm(v2)/layerThick)
   div=int(round(np.linalg.norm((v1+v2)/2.0)/layerThick))

   for lay in range(div):
      if (lay%2)==1:
         connector1_start=v1_start+(lay/div)*v1
         connector1=v2_start+(lay/div)*v2-connector1_start
         connector2_start=v1_start+((lay+1)/div)*v1
         connector2=v2_start+((lay+1)/div)*v2-connector2_start
      else:
         connector1_start=v2_start+(lay/div)*v2
         connector1=v1_start+(lay/div)*v1-connector1_start
         connector2_start=v2_start+((lay+1)/div)*v2
         connector2=v1_start+((lay+1)/div)*v1-connector2_start
      if abs(np.linalg.norm(connector1))>eps:
         connector1_start=projectOnVec(connector1_start,connector1,connector2_start)[0]
         connector1=-(connector1_start-projectOnVec(connector1_start,connector1,connector2_start+connector2)[0])
      increNum=int(np.linalg.norm(connector1)/increLen)
      lastp1=connector1_start
      lastp2=connector2_start
      for incre in range (increNum):
         p1=connector1_start+((incre+1)/increNum)*connector1
         p2=connector2_start+((incre+1)/increNum)*connector2
         #print([lastp1,lastp2,p1,p2])
         #ang1=np.arctan2(zFinder(p1)-zFinder(lastp1),abs(p1[0]-lastp1[0]))
         #ang2=np.arctan2(zFinder(p2)-zFinder(lastp2),abs(p2[0]-lastp2[0]))
         realDepth1=depth#/np.cos(ang1)
         realDepth2=depth#/np.cos(ang2)-
         cellCenter=(p1+p2+lastp1+lastp2)/4.0
         if not isInHole(cellCenter):
            pnt.write(str(round(lastp1[0],6))+' '+str(round(lastp1[1],6))+' '+str(round(zFinder(lastp1),6))+'\n')
            pnt.write(str(round(p1[0],6))+' '+str(round(p1[1],6))+' '+str(round(zFinder(p1),6))+'\n')
            pnt.write(str(round(p2[0],6))+' '+str(round(p2[1],6))+' '+str(round(zFinder(p2),6))+'\n')
            pnt.write(str(round(lastp2[0],6))+' '+str(round(lastp2[1],6))+' '+str(round(zFinder(lastp2),6))+'\n')
            pnt.write(str(round(lastp1[0],6))+' '+str(round(lastp1[1],6))+' '+str(round(zFinder(lastp1)-realDepth1,6))+'\n')
            pnt.write(str(round(p1[0],6))+' '+str(round(p1[1],6))+' '+str(round(zFinder(p1)-realDepth1,6))+'\n')
            pnt.write(str(round(p2[0],6))+' '+str(round(p2[1],6))+' '+str(round(zFinder(p2)-realDepth2,6))+'\n')
            pnt.write(str(round(lastp2[0],6))+' '+str(round(lastp2[1],6))+' '+str(round(zFinder(lastp2)-realDepth2,6))+'\n')
            elm.write("8 "+str(nCount)+' '+str(nCount+1)+' '+str(nCount+2)+' '+str(nCount+3)+' '+str(nCount+4)
                      +' '+str(nCount+5)+' '+str(nCount+6)+' '+str(nCount+7)+'\n')
            nCount=nCount+8
            elCount=elCount+1
         lastp1=p1
         lastp2=p2

############################################# OPTIMIZATION #########################################

def loopGB():
   #design parameters:
   global openRad
   openRad_start=0.2
   openRad_end=0.5
   
   global cenOpenX
   cenOpenX_start=3.0
   cenOpenX_end=1.6

   global cenOpenY
   cenOpenY_start=0.5
   cenOpenY_end=1.2

   global downAmp
   downAmp_start=0.0
   downAmp_end=0.1

   global downFreq
   downPeriod_start=1.5
   downPeriod_end=0.4

   global upAmp
   upAmp_start=0.0
   upAmp_end=0.08

   global upFreq
   upPeriod_start=1.0
   upPeriod_end=1.0

   global upPhase
   upPhase_start=0.0
   upPhase_end=np.pi/2.0

   global spotFreq
   spotPeriod_start=0.05
   spotPeriod_end=0.2

   global spotAmp
   spotAmp_start=0.0
   spotAmp_end=0.05

   #dependent variables
   global spotRad
   global spotEnd
   global spotX
   global spotY

   for ite in range(loopIters):
      openRad=openRad_start+(openRad_end-openRad_start)*ite/loopIters
      cenOpenX=cenOpenX_start+(cenOpenX_end-cenOpenX_start)*ite/loopIters
      cenOpenY=cenOpenY_start+(cenOpenY_end-cenOpenY_start)*ite/loopIters
      downAmp=downAmp_start+(downAmp_end-downAmp_start)*ite/loopIters
      downPeriod=downPeriod_start+(downPeriod_end-downPeriod_start)*ite/loopIters
      upAmp=upAmp_start+(upAmp_end-upAmp_start)*ite/loopIters
      upPeriod=upPeriod_start+(upPeriod_end-upPeriod_start)*ite/loopIters
      upPhase=upPhase_start+(upPhase_end-upPhase_start)*ite/loopIters
      spotPeriod=spotPeriod_start+(spotPeriod_end-spotPeriod_start)*ite/loopIters
      spotAmp=spotAmp_start+(spotAmp_end-spotAmp_start)*ite/loopIters
      downFreq=1.0/downPeriod
      upFreq=1.0/upPeriod
      spotFreq=1.0/spotPeriod
      spotRad=openRad
      spotEnd=openRad*2.0
      spotX=cenOpenX
      spotY=cenOpenY
      makeAwall(ite)
      
def loopDOE():
   #design parameters:
   global openRad
   openRad_min=0.2
   openRad_max=0.8
   
   global cenOpenX
   cenOpenX_min=0.5
   cenOpenX_max=4.0

   global cenOpenY
   cenOpenY_min=0.5
   cenOpenY_max=1.5

   global downAmp
   downAmp_min=0.0
   downAmp_max=0.25

   global downFreq
   downPeriod_min=0.2
   downPeriod_max=2.0

   global upAmp
   upAmp_min=0.0
   upAmp_max=0.25

   global upFreq
   upPeriod_min=0.2
   upPeriod_max=2.0

   global upPhase
   upPhase_min=0.0
   upPhase_max=np.pi/2.0

   global spotFreq
   spotPeriod_min=0.1
   spotPeriod_max=0.25

   global spotAmp
   spotAmp_min=0.0
   spotAmp_max=0.15

   #dependent variables
   global spotRad
   global spotEnd
   global spotX
   global spotY

   for ite in range(loopIters):
      openRad=rand.uniform(openRad_min, openRad_max)
      cenOpenX=rand.uniform(cenOpenX_min, cenOpenX_max)
      cenOpenY=rand.uniform(cenOpenY_min, cenOpenY_max)
      downAmp=rand.uniform(downAmp_min, downAmp_max)
      downPeriod=rand.uniform(downPeriod_min, downPeriod_max)
      upAmp=rand.uniform(upAmp_min, upAmp_max)
      upPeriod=rand.uniform(upPeriod_min, upPeriod_max)
      upPhase=rand.uniform(upPhase_min, upPhase_max)
      spotPeriod=rand.uniform(spotPeriod_min, spotPeriod_max)
      spotAmp=rand.uniform(spotAmp_min, spotAmp_max)
      downFreq=1.0/downPeriod
      upFreq=1.0/upPeriod
      spotFreq=1.0/spotPeriod
      spotRad=openRad
      spotEnd=openRad*2.0
      spotX=cenOpenX
      spotY=cenOpenY
      makeAwall(ite)

############################################# OUTPUT WRITING #########################################

def plotWall():
   fig = pyplot.figure()
   ax = Axes3D(fig)
   ax.scatter(x,y,z)
   ax.set_xlim3d(0,max(wallWidth,wallHeight))
   ax.set_ylim3d(0,max(wallWidth,wallHeight))
   pyplot.show()
   
def writeCSV():
   with open('wallPoints.csv', mode='w') as pointsFile:
      pointWriter = csv.writer(pointsFile, delimiter=',')
      for i in range(xSize):
         for j in range(ySize):
            pointWriter.writerow([x[i][j],y[i][j],z[i][j]])

def writeVTK(ite):
   if withHole:
      findMissingEls()
   numMissing=0
   for el in elMissing:
      if el==1:
          numMissing=numMissing+1
   counter=0
   vtk = open("wallPoints.vtk", "w")
   vtk.write("# vtk DataFile Version 1.0\nUnstructured Grid Example\nASCII\n\nDATASET UNSTRUCTURED_GRID\nPOINTS ")
   vtk.write(str(xSize*ySize)+" float\n")
   for i in range(xSize):
      for j in range(ySize):
         vtk.write(str(round(x[i][j],6))+' '+str(round(y[i][j],6))+' '+str(round(z[i][j],6))+'\n')
         counter=counter+1
   vtk.write("\nCELLS "+str((xSize-1)*(ySize-1)-numMissing)+' '+str(5*((xSize-1)*(ySize-1)-numMissing))+'\n')
   counter=0
   for i in range(xSize-1):
      for j in range(ySize-1):
         if elMissing[counter]==0:
            vtk.write("4 "+str(i*ySize+j)+' '+str(i*ySize+j+1)+' '+str((i+1)*ySize+j+1)+' '+str((i+1)*ySize+j)+'\n')
         counter=counter+1
   counter=0
   vtk.write("\nCELL_TYPES "+str((xSize-1)*(ySize-1)-numMissing)+'\n')
   for i in range(xSize-1):
      for j in range(ySize-1):
         if elMissing[counter]==0:
             vtk.write("9\n")
         counter=counter+1

   if withCurv:
      counter=0
      vtk.write("\nPOINT_DATA "+str(xSize*ySize)+"\nSCALARS GaussCurv float\nLOOKUP_TABLE default\n")
      for i in range(xSize):
         for j in range(ySize):
            vtk.write(str(round(curv[i][j],6))+'\n')
            counter=counter+1

   if withOverAngle:      
      counter=0
      vtk.write("\nSCALARS OverHang float\nLOOKUP_TABLE default\n")
      for i in range(xSize):
         for j in range(ySize):
            vtk.write(str(round(overAngle[i][j],6))+'\n')
            counter=counter+1

   if withGate:      
      counter=0
      vtk.write("\nVECTORS disp float\n")
      for i in range(xSize):
         for j in range(ySize):
            vtk.write(str(round(dispX[i][j],6))+' '+str(round(dispY[i][j],6))+" 0.0 \n")
            counter=counter+1   

def writeAMvtk():
   dirName="AMpathFrames/"
   if not(os.path.exists(dirName)):
      os.mkdir(dirName)
   for i in range (frames):
      writeVTKframe(i*int(elCount/frames),i,dirName)
   if elCount>(frames*int(elCount/frames)):
      writeVTKframe(elCount,frames,dirName)


def writeVTKframe(psTime,frame,dirName):

   fileName="AMpath"+str(frame)+".vtk"
   vtk = open(dirName+fileName, "w")
   pnt = open("AMpath.pnt", "r")
   elm = open("AMpath.elm", "r")

   vtk.write("# vtk DataFile Version 1.0\nUnstructured Grid Example\nASCII\n\nDATASET UNSTRUCTURED_GRID\nPOINTS ")
   vtk.write(str(psTime*8)+" float\n")
   pntLines = pnt.readlines()
   for line in range(psTime*8): 
      vtk.write(pntLines[line])
   vtk.write("\nCELLS "+str(psTime)+' '+str(psTime*9)+'\n')
   elmLines = elm.readlines()
   for line in range(psTime): 
      vtk.write(elmLines[line].strip()+"\n")
   vtk.write("\nCELL_TYPES "+str(psTime)+"\n")
   for i in range (psTime):
      vtk.write("12\n")

   if withDrieness:
      timePerCell=nozzleSpeed*increLen
      vtk.write("\nCELL_DATA "+str(psTime)+"\nSCALARS age float\nLOOKUP_TABLE default\n")
      for line in range(psTime): 
         vtk.write(str(round(timePerCell*(psTime-line),6))+"\n")
      
   elm.close()
   pnt.close()
   vtk.close()
   
main()
   



