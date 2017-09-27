"""
Example 2: Platooning
    Here we find optimal communication graphs for a platoon
"""

import random
import numpy as np
from ana_invariance import *

s=system()
Np=6 # Number of platoons
s.n=2*Np # number of variables
s.m=Np # number of controls
s.K=8 # Design variable
tau=0.1 # coefficient of disturbances for x


varepsilon=0.1836
beta=varepsilon

for i in range(1,s.n+1):
    for j in range(1,s.n+1):
        s.A[i,j]=int(i==j)
        
for i in range(2,Np+1):
    s.A[2*i-1,2*i]=-1
    s.A[2*i-1,2*i-2]=1

s.A[1,2]=-1   

for i in range(1,s.n+1):
    for j in range(1,s.m+1):
        s.B[i,j]=0

for i in range(1,Np+1):
    s.B[2*i,i]=1


print "\nA=\n",matricize(s.A)
print "\nB=\n",matricize(s.B)


for i in range(1,Np+2):
    for j in range(1,s.n+1):
        s.Hx[i,j]=0


for i in range(1,Np+1):
    s.Hx[i,2*i-1]=-1
    
for i in range(1,Np+1):
    s.Hx[Np+1,2*i-1]=1

for i in range(1,Np+1):
    s.hx[i]=0.5
s.hx[Np+1]=Np/2.0
    
print "\nHx=\n",matricize(s.Hx)
print "\nhx=\n",s.hx


# Constructing W:
for i in range(1,4*Np+Np*(Np-1)+1):
    for j in range(1,s.n+1):
        s.Hw[i,j]=0

# The part that is easy: squares        
for i in range(1,s.n+1):
    s.Hw[2*i-1,i]=-1
    s.Hw[2*i,i]=1

for i in range(1,Np+1):
    s.hw[4*i-3]=varepsilon*tau
    s.hw[4*i-2]=varepsilon*tau
    s.hw[4*i-1]=varepsilon+beta
    s.hw[4*i]=varepsilon+beta

r=4*Np+1    
for i in range(1,Np+1):
    for j in range(i+1,Np+1):
        s.Hw[r,2*i]=1
        s.Hw[r,2*j]=-1
        s.Hw[r+1,2*i]=-1
        s.Hw[r+1,2*j]=1
        r+=2

for i in range(4*Np+1,4*Np+Np*(Np-1)+1):
    s.hw[i]=2*varepsilon


print "\nHw=\n",matricize(s.Hw)
print "hw=",s.hw

for i in range(1,2*s.m+1):
    for j in range(1,s.m+1):
        s.Hu[i,j]=0


for i in range(1,s.m+1):
    s.Hu[2*i-1,i]=1
    s.Hu[2*i,i]=-1
    

for i in range(1,2*s.m+1):
    s.hu[i]=1

print "\nHu=\n",matricize(s.Hu)
print "\nhu=\n",s.hu
            
s.N=Np+1
for i in range(1,Np+1):
    s.subsystem_X[2*i-1]=i+1
    s.subsystem_X[2*i]=i+1
    s.subsystem_U[i]=i+1

print "s.X:",s.subsystem_X
print "s.U:",s.subsystem_U

for i in range(1,Np+2):
    for j in range(1,Np+2):
        s.c[i,j]=abs(i-j)**2

# s.RCI()
s.GraphDesign() # Computes a graph

print "\nCommunications=\n",matricize(s.C)


print "nX=", s.nX
print "nU=", s.nU
print "nW", s.nW

# raise "stop!"


if False: # Raise it as false
    f=open("RCI.txt","w")
    print "RCI set Projection"
    x={}
    for i in range(1,s.n+1):
        x[i]=0
    N=200
    for p in range(0,N):
        print "row ",p
        for q in range(0,N):
            x1=-1+p/float(N)*2
            x2=-1+q/float(N)*2
            x[1]=x1
            x[2]=x2
            f.write("%d "%int(s.is_RCI(x)))
        f.write("\n")
    f.close()           
        
    
    
    
    

if False:
    print '-'*80
    print "\t\t\t\tSimulation:"
    print '-'*80
    # A test with new idea
    for t in range(0,s.K+1):
        for i in range(1,s.n+1):
            s.xH[t,i]=0
    for t in range(0,s.K):
        for i in range(1,s.m+1):
            s.uH[t,i]=0        
    s.t=s.K
    f=open("platoon_state.txt","w")
    for t in range(s.K,60):
        w={}
        w[1]=0.0*random.random()**0.2*(-1)**(random.random()>0.5)
        for i in range(2,s.n+1):
            if i%2!=0:
                w[i]=varepsilon*random.random()**0.2*(-1)**(random.random()>0.5) 
            else:
                w[i]=0
        print t,"x[1]:",s.xH[t,1],"x[2]:",s.xH[t,2],"x[3]:",s.xH[t,3],"x[2]:",s.xH[t,4],"x[5]:",s.xH[t,5]
        s.compute_control_delay()
        print t,"u[1]:",s.u[1],"u[2]:",s.u[2],"u[3]:",s.u[3]
        print "**********"
        for i in range(1,s.n+1):
            f.write("%0.2f " %s.xH[t,i])
        f.write("\n")
        s.apply_control_delay(w)
        print "is state in RCI set?",s.is_RCI(s.x)
    f.close()
else:
    print "I'm not running the simulation"
    



    