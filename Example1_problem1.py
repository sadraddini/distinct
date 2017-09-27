"""
Example 1
"""

import random
import numpy as np
from ana_invariance import *

s=system()
s.n=10 # number of variables
s.m=5 # number of controls
s.K=6 # Design variable
c=0 # 0 for directed, 1 for undirected

eps=0.01
eta=0.1

for i in range(1,s.n+1):
    for j in range(1,s.n+1):
        s.A[i,j]=eps*(-1)**(i+j)
        s.A[i,j]+=int(i==j)
        

s.A[1,2]=1
s.A[3,4]=1
s.A[5,6]=1
s.A[7,8]=1
s.A[9,10]=1

for i in range(1,11):
    for j in range(1,6):
        s.B[i,j]=eps*(-1)**(i+2*j)

for i in range(1,s.m+1):
    s.B[2*i-1,i]=0
    s.B[2*i,i]=1


for i in range(1,21):
    for j in range(1,11):
        s.Hw[i,j]=0


for i in range(1,s.n+1):
    s.Hw[2*i-1,i]=1
    s.Hw[2*i,i]=-1


for i in range(1,2*s.n+1):
    s.hw[i]=eta


s.Hx=s.Hw

for i in range(1,2*s.n+1):
    s.hx[i]=1

for i in range(1,11):
    for j in range(1,6):
        s.Hu[i,j]=0


for i in range(1,s.m+1):
    s.Hu[2*i-1,i]=1
    s.Hu[2*i,i]=-1
    

for i in range(1,2*s.m+1):
    s.hu[i]=2


print "\nA=\n",matricize(s.A)
print "\nB=\n",matricize(s.B)
print "\nHx=\n",matricize(s.Hx)
print "\nHu=\n",matricize(s.Hu)
print "\nHw=\n",matricize(s.Hw)
				
# s.RCI() # Computes a set
# raise "stop!"
s.N=5
for i in range(1,s.N+1):
    s.subsystem_X[2*i-1]=i
    s.subsystem_X[2*i]=i
    s.subsystem_U[i]=i

for i in range(1,s.N+1):
    for j in range(1,s.N+1):
        s.G[i,j,1]=int(i==j)
        
s.G[1,2,1]=1
s.G[2,1,1]=c

s.G[2,3,1]=1
s.G[3,2,1]=c

s.G[3,4,1]=1
s.G[4,3,1]=c

s.G[4,5,1]=1
s.G[5,4,1]=c

s.G[5,1,1]=1
s.G[1,5,1]=c
    
s.graph_powers()
for k in range(1,s.K+2):
    print "G",k,":["
    for i in range(1,s.N+1):
        for j in range(1,s.N+1):
            print s.G[i,j,k],
        print ""
    print "]"

# raise "stop"
s.Margin_Design() # Computes a distributed-based set

print "nX=", s.nX
print "nU=", s.nU
print "nW", s.nW



if True: # Raise it as false
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
        
    
    
    
    

if True:
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
    f=open("state.txt","w")
    for t in range(s.K,60):
        w={}
        for i in range(1,s.n+1):
            w[i]=eta*random.random()**0.2*(-1)**(random.random()>0.5) 
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
    



	