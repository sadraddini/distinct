""" 
Version 1.0

    Synthesizing Distributed Set-Invariance Policies for Linear 
    Discrete-Time Systems Subject to Additive Disturbances

Based on paper:
    Sadra Sadraddini, and Calin Belta, 
    "Distributed Set-Invariance Control for Interconnected Linear Systems"
    Submitted to American Control Conference, 2018

Developed by:
    Sadra Sadraddini
    Boston University, Boston, MA
    sadra at bu dot edu
    September 2017

Notes:
....... Usage with proper citation is allowed
....... The code is still under development
"""

from gurobipy import *
import numpy as np
import random

class system:
    def __init__(self):
        self.A={} #  A matrix, enter as a dictionary
        self.B={} #  B matrix, enter as a dictionary
        
        self.Hx={} #  X={x| Hx x\le hx}
        self.hx={}
        self.Hu={} #  U={u| Hu u\le hu}
        self.hu={}
        self.Hw={} #  W={w| Hw w\le hw}
        self.hw={}           
        
        self.alpha=0.0 # Contraction Factor
        
        self.n=0 # number of variables
        self.m=0 # number of controls
        self.K=1 # Design variable, degree
        self.nX=1 # rows of Hx 
        self.nU=1 # rows of Hu
        self.nW=1 # rows of Hw

        self.AA={} # AA[k,i,j] is i,j entry of A^k
        self.HxAB={} # HxAB[k,i,j] is i,j entry of Hx*A^k*B
        self.HwAB={} # HwAB[k,i,j] is i,j entry of Hw*A^k*B
        self.theta={}
        self.D={}
        
        # Dynamic Variables:
        self.x={}
        self.u={}
        self.w={}
        
        # History Variables: xH[t,i]: x[i] at time t
        self.xH={}
        self.uH={}
        self.wH={}
        self.t=0 # Time!
        
        # How subsystems are distinguished
        self.N=1 # Number of subsystems
        self.subsystem_X={} # subsystem_X[i]=s: i belongs to s
        self.subsystem_U={} # subsystem_U[i]=s: i belongs to s
        self.G={} # Graph Variables
        self.C={} # Communication Structure
        self.c={} # The costs!
        
    def find_rows(self):
        self.nX=max(self.hx.keys())
        self.nU=max(self.hu.keys())
        self.nW=max(self.hw.keys())
        
    def compute_HxAB(self):
        for k in range(0,self.K):
            for i in range(1,self.nX+1):
                for j in range(1,self.m+1):
                    self.HxAB[k,i,j]=0
                    for q in range(1,self.n+1):
                        for p in range(1,self.n+1):
                            self.HxAB[k,i,j]+=self.Hx[i,q]*self.AA[k,q,p]*self.B[p,j]
    def compute_HwAB(self):
        for k in range(0,self.K):
            for i in range(1,self.nW+1):
                for j in range(1,self.m+1):
                    self.HwAB[k,i,j]=0
                    for q in range(1,self.n+1):
                        for p in range(1,self.n+1):
                            self.HwAB[k,i,j]+=self.Hw[i,q]*self.AA[k,q,p]*self.B[p,j]
    def compute_AA(self):
        for i in range(1,self.n+1):
            for j in range(1,self.n+1):
                self.AA[0,i,j]=int(i==j)
                self.AA[1,i,j]=self.A[i,j]
        for k in range(0,self.K):
            for i in range(1,self.n+1):
                for j in range(1,self.n+1):
                    self.AA[k+1,i,j]=0
                    for p in range(1,self.n+1):
                        self.AA[k+1,i,j]+=self.A[i,p]*self.AA[k,p,j]
                
    def compute_D(self):
        for k in range(0,self.K):
            for i in range(1,self.n+1):
                for j in range(1,self.n+1):
                    self.D[k,i,j]=self.AA[self.K-k-1,i,j]
                    for p in range(1,self.n+1):
                        for q in range(1,self.m+1):
                            for eta in range(0,self.K-k-1):
                                self.D[k,i,j]+=self.AA[self.K-k-2-eta,i,p]*self.B[p,q]*self.theta[eta,q,j]
                    if abs(self.D[k,i,j])<10**(-13):
                        self.D[k,i,j]=0
 
    def RCI(self):
        M={}
        GAM={}
        PI={}
        LAMBDA={}
        xbar={}
        ubar={}
        model=Model("RCIS")
        self.find_rows()
        self.compute_AA()
        self.compute_HxAB()
        self.compute_HwAB()
        # These are new!
        self.beta=model.addVar(lb=0,ub=1)
        self.gamma=model.addVar(lb=0,ub=1)
        self.rho=model.addVar(lb=0,ub=GRB.INFINITY,obj=1)
        
        
        for i in range(1,self.n+1):
            xbar[i]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
        for i in range(1,self.m+1):
            ubar[i]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
            
        for k in range(0,self.K):
            for i in range(1,self.m+1):
                for j in range(1,self.n+1):
                    self.theta[k,i,j]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
        for i in range(1,self.nW+1):
            for j in range(1,self.nW+1):
                LAMBDA[i,j]=model.addVar(lb=0)
                
        for k in range(0,self.K):
            for i in range(1,self.nX+1):
                for j in range(1,self.nW+1):
                    GAM[k,i,j]=model.addVar(lb=0,ub=GRB.INFINITY)
            for i in range(1,self.nU+1):
                for j in range(1,self.nW+1):
                    PI[k,i,j]=model.addVar(lb=0,ub=GRB.INFINITY)        
            
        model.update()
        
        # A xbar + B u bar= xbar
        for i in range(1,self.n+1):
            s=LinExpr()
            for j in range(1,self.n+1):
                s.addTerms(self.A[i,j], xbar[j])
            for j in range(1,self.m+1):
                s.addTerms(self.B[i,j], ubar[j])
            model.addConstr(s == xbar[i])
        
        for i in range(1,self.n+1):
            model.addConstr(xbar[i]==0)
        for j in range(1,self.m+1):
            model.addConstr(ubar[j]==0)
        
        # Equation 4.9-1
        # Lambda * hw <= alpha * hw
        for i in range(1,self.nW+1):
            s=LinExpr()
            for j in range(1,self.nW+1):
                s.addTerms(self.hw[j],LAMBDA[i,j])
            model.addConstr(s <= self.alpha*self.hw[i])
        
        # Equation 4.9-2
        for i in range(1,self.nW+1):
            for j in range(1,self.n+1):
                leftHand=LinExpr()
                for p in range(1,self.nW+1):
                    leftHand.addTerms(self.Hw[p,j],LAMBDA[i,p])
                # Now the heavy part
                rightHand=LinExpr()
                for eta in range(0,self.K):
                    for p in range(1,self.m+1): 
                        rightHand.addTerms(self.HwAB[self.K-1-eta,i,p],self.theta[eta,p,j])
                righthandConstant=0
                for p in range(1,self.n+1):
                    righthandConstant+=self.Hw[i,p]*self.AA[self.K,p,j]
                model.addConstr( leftHand == rightHand + righthandConstant)
        
        # Equation 4.10a)
        for i in range(1,self.nX+1):
            s=LinExpr()
            Hx=LinExpr()
            for j in range(1,self.n+1):
                Hx.addTerms(-self.Hx[i,j],xbar[j])
            for k in range(0,self.K):
                for j in range(1,self.nW+1):
                    s.addTerms(self.hw[j],GAM[k,i,j])
            model.addConstr(s <= (1-self.alpha)*self.beta*self.hx[i] + (1-self.alpha)*Hx )
          
        # Equation 4.10b)
        for k in range(0,self.K):
            for i in range(1,self.nX+1):
                for j in range(1,self.n+1):
                    leftHand=LinExpr()
                    for p in range(1,self.nW+1):
                        leftHand.addTerms(self.Hw[p,j],GAM[k,i,p])
                    # Now the heavy part
                    rightHand=LinExpr()
                    for eta in range(0,k):
                        for p in range(1,self.m+1): 
                            rightHand.addTerms(self.HxAB[k-1-eta,i,p],self.theta[eta,p,j])
                    righthandConstant=0
                    for p in range(1,self.n+1):
                        righthandConstant+=self.Hx[i,p]*self.AA[k,p,j]
                    model.addConstr( leftHand == rightHand + righthandConstant)
        
        # Equation 4.11a)
        for i in range(1,self.nU+1):
            s=LinExpr()
            Pu=LinExpr()
            for j in range(1,self.m+1):
                Pu.addTerms(-self.Hu[i,j],ubar[j])
            for k in range(0,self.K):
                for j in range(1,self.nW+1):
                    s.addTerms(self.hw[j],PI[k,i,j])
            model.addConstr(s <= (1-self.alpha)*self.gamma*self.hu[i]+ (1-self.alpha)*Pu )
        
        # Equation 4.11b)
        for k in range(0,self.K):
            for i in range(1,self.nU+1):
                for j in range(1,self.n+1):
                    leftHand=LinExpr()
                    for p in range(1,self.nW+1):
                        leftHand.addTerms(self.Hw[p,j],PI[k,i,p])
                    # Now the heavy part
                    rightHand=LinExpr()
                    for p in range(1,self.m+1): 
                        rightHand.addTerms(self.Hu[i,p],self.theta[k,p,j])
                    model.addConstr( leftHand == rightHand)            
        
        model.addConstr( self.rho >= self.beta)
        model.addConstr( self.rho >= self.gamma)    
        
        model.optimize()
        if model.Status==3:
            print "RCI is empty"
            return False
        else:
            print "\n\tGreat! RCI set exists and was computed!\n"
            for k in range(0,self.K):
                for i in range(1,self.m+1):
                    for j in range(1,self.n+1):
                        self.theta[k,i,j]=self.theta[k,i,j].X
            print "beta was", self.beta.X
            print "gamma was", self.gamma.X
            self.compute_D()
            return True
                    
    
    
    
    
    
    
    def graph_powers(self):
        """ 
            Here I design powers of G
            I assume values for G are given, otherwise error will be returned
        """
        for k in range(1,self.K+1):
            for I in range(1,self.N+1):
                for J in range(1,self.N+1):
                    self.G[I,J,k+1]=self.G[I,J,k]
                    for p in range(1,self.N+1):
                        self.G[I,J,k+1]+=self.G[I,p,k]*self.G[p,J,1]
                    self.G[I,J,k+1]=min(self.G[I,J,k+1],1)
                        
                    
            
        
    
    
    
    
    def Margin_Design(self):
        """
            Problem 1 in the paper
            Here I want to maximize the margin of correctness
            The graph is assumed to be given
            The powers should have been computed before hand!
        """
        GAM={}
        PI={}
        LAMBDA={}
        model=Model("RCIS")
        self.find_rows()
        self.compute_AA()
        self.compute_HwAB()
        self.compute_HxAB()
        self.graph_powers()
        # These are implemented!
        self.beta=model.addVar(lb=0,ub=1)
        self.gamma=model.addVar(lb=0,ub=1)
        self.rho=model.addVar(lb=0,ub=1,obj=1)               
            
        for k in range(0,self.K):
            for i in range(1,self.m+1):
                for j in range(1,self.n+1):
                    self.theta[k,i,j]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
        for i in range(1,self.nW+1):
            for j in range(1,self.nW+1):
                LAMBDA[i,j]=model.addVar(lb=0)
                
        for k in range(0,self.K):
            for i in range(1,self.nX+1):
                for j in range(1,self.nW+1):
                    GAM[k,i,j]=model.addVar(lb=0,ub=GRB.INFINITY)
            for i in range(1,self.nU+1):
                for j in range(1,self.nW+1):
                    PI[k,i,j]=model.addVar(lb=0,ub=GRB.INFINITY)        
            
        model.update()
        
        # Equation 4.9-1
        # Lambda * g <= alpha * g
        for i in range(1,self.nW+1):
            s=LinExpr()
            for j in range(1,self.nW+1):
                s.addTerms(self.hw[j],LAMBDA[i,j])
            model.addConstr(s <= self.alpha*self.hw[i])
        
        # Equation 4.9-2
        for i in range(1,self.nW+1):
            for j in range(1,self.n+1):
                leftHand=LinExpr()
                for p in range(1,self.nW+1):
                    leftHand.addTerms(self.Hw[p,j],LAMBDA[i,p])
                # Now the heavy part
                rightHand=LinExpr()
                for eta in range(0,self.K):
                    for p in range(1,self.m+1): 
                        rightHand.addTerms(self.HwAB[self.K-1-eta,i,p],self.theta[eta,p,j])
                righthandConstant=0
                for p in range(1,self.n+1):
                    righthandConstant+=self.Hw[i,p]*self.AA[self.K,p,j]
                model.addConstr( leftHand == rightHand + righthandConstant)
        
        # Equation 4.10a)
        for i in range(1,self.nX+1):
            s=LinExpr()
            for k in range(0,self.K):
                for j in range(1,self.nW+1):
                    s.addTerms(self.hw[j],GAM[k,i,j])
            model.addConstr(s <= (1-self.alpha)*self.beta*self.hx[i])
          
        # Equation 4.10b)
        for k in range(0,self.K):
            for i in range(1,self.nX+1):
                for j in range(1,self.n+1):
                    leftHand=LinExpr()
                    for p in range(1,self.nW+1):
                        leftHand.addTerms(self.Hw[p,j],GAM[k,i,p])
                    # Now the heavy part
                    rightHand=LinExpr()
                    for eta in range(0,k):
                        for p in range(1,self.m+1): 
                            rightHand.addTerms(self.HxAB[k-1-eta,i,p],self.theta[eta,p,j])
                    righthandConstant=0
                    for p in range(1,self.n+1):
                        righthandConstant+=self.Hx[i,p]*self.AA[k,p,j]
                    model.addConstr( leftHand == rightHand + righthandConstant)
        
        # Equation 4.11a)
        for i in range(1,self.nU+1):
            s=LinExpr()
            for k in range(0,self.K):
                for j in range(1,self.nW+1):
                    s.addTerms(self.hw[j],PI[k,i,j])
            model.addConstr(s <= (1-self.alpha)*self.gamma*self.hu[i] )
        
        # Equation 4.11b)
        for k in range(0,self.K):
            for i in range(1,self.nU+1):
                for j in range(1,self.n+1):
                    leftHand=LinExpr()
                    for p in range(1,self.nW+1):
                        leftHand.addTerms(self.Hw[p,j],PI[k,i,p])
                    # Now the heavy part
                    rightHand=LinExpr()
                    for p in range(1,self.m+1): 
                        rightHand.addTerms(self.Hu[i,p],self.theta[k,p,j])
                    model.addConstr( leftHand == rightHand)                
        # Control Structure
        """
            \theta_k B \in G^{k+1}, k=0 to K-1
        """
        for k in range(0,self.K):
            for i in range(1,self.m+1):
                for j in range(1,self.m+1):
                    I=self.subsystem_U[i]
                    J=self.subsystem_U[j]
                    s=LinExpr()
                    for p in range(1,self.n+1):
                        s.add(self.theta[k,i,p]*self.B[p,j])
                    if self.G[J,I,k+1]==0:
                        model.addConstr(s==0)
        """
         State Structure
         1) M0 \in G1
         2) M(k)-M(k-1)A \in G(k+1)
         3) M(K-1) A \in G(K+1)
        """
        # 1: M0 \in G1
        for i in range(1,self.m+1):
            for j in range(1,self.n+1):
                I=self.subsystem_U[i]
                J=self.subsystem_X[j]
                if self.G[J,I,1]==0:
                    model.addConstr( self.theta[0,i,j] == 0)
        
        # 2: M(k)-M(k-1)A \in G(k+1)
        for k in range(1,self.K):
            for i in range(1,self.m+1):
                for j in range(1,self.n+1):
                    I=self.subsystem_U[i]
                    J=self.subsystem_X[j]
                    s=LinExpr()
                    s.add(self.theta[k,i,j])
                    for p in range(1,self.n+1):
                        s.add(-self.theta[k-1,i,p]*self.A[p,j]) 
                    if self.G[J,I,k+1]==0:
                        model.addConstr( s == 0)

        # 3: M(k-1)*A \in G(K+1)
        for i in range(1,self.m+1):
            for j in range(1,self.n+1):
                I=self.subsystem_U[i]
                J=self.subsystem_X[j]
                s=LinExpr()
                for p in range(1,self.n+1):
                    s.add(self.theta[self.K-1,i,p]*self.A[p,j]) 
                if self.G[J,I,self.K+1]==1:
                    model.addConstr( s==0)
        
        model.addConstr( self.rho >= self.beta)
        model.addConstr( self.rho >= self.gamma)
        
        model.optimize()
        if model.Status==3:
            print "RCI is empty"
            return False
        else:
            print "\n\tGreat! RCI set exists and was computed!\n"
            for k in range(0,self.K):
                for i in range(1,self.m+1):
                    for j in range(1,self.n+1):
                        self.theta[k,i,j]=self.theta[k,i,j].X
            self.compute_D()
            print "The margin of correctness was", 1.00-self.rho.X
            return True    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    def GraphDesign(self):
        """
            Problem 2 in the paper
            Here I want to minimize total communication cost
        """
        GAM={}
        PI={}
        LAMBDA={}
        model=Model("RCIS")
        self.find_rows()
        self.compute_AA()
        self.compute_HwAB()
        self.compute_HxAB()
        
        # These are NOT implemented!
        self.beta=model.addVar(lb=0,ub=1)
        self.gamma=model.addVar(lb=0,ub=1)
        self.rho=model.addVar(lb=0,ub=1,obj=0)
        
        # Graph representations
        for k in range(1,self.K+2):
            for i in range(1,self.N+1):
                for j in range(1,self.N+1):
                    if i==j:
                        self.G[i,j,k]=1
                    elif k==1:
                        self.G[i,j,1]=model.addVar(vtype=GRB.BINARY,obj=self.c[i,j])
                    else:
                        self.G[i,j,k]=model.addVar(lb=0,ub=1)
                    for p in range(1,self.N+1):
                        self.G["temp",p,i,j,k]=model.addVar(lb=0,ub=1)               
            
        for k in range(0,self.K):
            for i in range(1,self.m+1):
                for j in range(1,self.n+1):
                    self.theta[k,i,j]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
        for i in range(1,self.nW+1):
            for j in range(1,self.nW+1):
                LAMBDA[i,j]=model.addVar(lb=0)
                
        for k in range(0,self.K):
            for i in range(1,self.nX+1):
                for j in range(1,self.nW+1):
                    GAM[k,i,j]=model.addVar(lb=0,ub=GRB.INFINITY)
            for i in range(1,self.nU+1):
                for j in range(1,self.nW+1):
                    PI[k,i,j]=model.addVar(lb=0,ub=GRB.INFINITY)        
            
        model.update()
        
        # Equation 4.9-1
        # Lambda * g <= alpha * g
        for i in range(1,self.nW+1):
            s=LinExpr()
            for j in range(1,self.nW+1):
                s.addTerms(self.hw[j],LAMBDA[i,j])
            model.addConstr(s <= self.alpha*self.hw[i])
        
        # Equation 4.9-2
        for i in range(1,self.nW+1):
            for j in range(1,self.n+1):
                leftHand=LinExpr()
                for p in range(1,self.nW+1):
                    leftHand.addTerms(self.Hw[p,j],LAMBDA[i,p])
                # Now the heavy part
                rightHand=LinExpr()
                for eta in range(0,self.K):
                    for p in range(1,self.m+1): 
                        rightHand.addTerms(self.HwAB[self.K-1-eta,i,p],self.theta[eta,p,j])
                righthandConstant=0
                for p in range(1,self.n+1):
                    righthandConstant+=self.Hw[i,p]*self.AA[self.K,p,j]
                model.addConstr( leftHand == rightHand + righthandConstant)
        
        # Equation 4.10a)
        for i in range(1,self.nX+1):
            s=LinExpr()
            for k in range(0,self.K):
                for j in range(1,self.nW+1):
                    s.addTerms(self.hw[j],GAM[k,i,j])
            model.addConstr(s <= (1-self.alpha)*self.beta*self.hx[i])
          
        # Equation 4.10b)
        for k in range(0,self.K):
            for i in range(1,self.nX+1):
                for j in range(1,self.n+1):
                    leftHand=LinExpr()
                    for p in range(1,self.nW+1):
                        leftHand.addTerms(self.Hw[p,j],GAM[k,i,p])
                    # Now the heavy part
                    rightHand=LinExpr()
                    for eta in range(0,k):
                        for p in range(1,self.m+1): 
                            rightHand.addTerms(self.HxAB[k-1-eta,i,p],self.theta[eta,p,j])
                    righthandConstant=0
                    for p in range(1,self.n+1):
                        righthandConstant+=self.Hx[i,p]*self.AA[k,p,j]
                    model.addConstr( leftHand == rightHand + righthandConstant)
        
        # Equation 4.11a)
        for i in range(1,self.nU+1):
            s=LinExpr()
            for k in range(0,self.K):
                for j in range(1,self.nW+1):
                    s.addTerms(self.hw[j],PI[k,i,j])
            model.addConstr(s <= (1-self.alpha)*self.gamma*self.hu[i] )
        
        # Equation 4.11b)
        for k in range(0,self.K):
            for i in range(1,self.nU+1):
                for j in range(1,self.n+1):
                    leftHand=LinExpr()
                    for p in range(1,self.nW+1):
                        leftHand.addTerms(self.Hw[p,j],PI[k,i,p])
                    # Now the heavy part
                    rightHand=LinExpr()
                    for p in range(1,self.m+1): 
                        rightHand.addTerms(self.Hu[i,p],self.theta[k,p,j])
                    model.addConstr( leftHand == rightHand)                
        
        
        # Graph Equations:
        for k in range(1,self.K):
            for i in range(1,self.N+1):
                for j in range(1,self.N+1):
                    s=LinExpr()
                    s.add(self.G[i,j,k])
                    if i!=j:
                        model.addConstr( self.G[i,j,k+1]>=self.G[i,j,k])
                    for p in range(1,self.N+1):
                        model.addConstr( self.G["temp",p,i,j,k] <= self.G[i,p,k] )
                        model.addConstr( self.G["temp",p,i,j,k] <= self.G[p,j,k] )
                        model.addConstr( self.G["temp",p,i,j,k] >= self.G[i,p,k]+self.G[p,j,k]-1)
                        s.add(self.G["temp",p,i,j,k])
                        model.addConstr( self.G[i,j,k+1] >= self.G["temp",p,i,j,k] )
                    model.addConstr( self.G[i,j,k+1] <= s)
        
        bigM=10
        # Control Structure
        """
            \theta B
        """
        for k in range(0,self.K):
            for i in range(1,self.m+1):
                for j in range(1,self.m+1):
                    I=self.subsystem_U[i]
                    J=self.subsystem_U[j]
                    s=LinExpr()
                    for p in range(1,self.n+1):
                        s.add(self.theta[k,i,p]*self.B[p,j])
                    model.addConstr( s <= bigM * self.G[J,I,k+1])
                    model.addConstr( s >= -bigM * self.G[J,I,k+1])
        """
         State Structure
         1) M0 \in G1
         2) M(k)-M(k-1)A \in G(k+1)
         3) M(K-1) A \in G(K+1)
        """
        # 1: M0 \in G1
        for i in range(1,self.m+1):
            for j in range(1,self.n+1):
                I=self.subsystem_U[i]
                J=self.subsystem_X[j]
                model.addConstr( self.theta[0,i,j] <= bigM * self.G[J,I,1])
                model.addConstr( self.theta[0,i,j] >= -bigM * self.G[J,I,1])
        
        # 2: M(k)-M(k-1)A \in G(k+1)
        for k in range(1,self.K):
            for i in range(1,self.m+1):
                for j in range(1,self.n+1):
                    I=self.subsystem_U[i]
                    J=self.subsystem_X[j]
                    s=LinExpr()
                    s.add(self.theta[k,i,j])
                    for p in range(1,self.n+1):
                        s.add(-self.theta[k-1,i,p]*self.A[p,j]) 
                    model.addConstr( s <= bigM * self.G[J,I,k+1])
                    model.addConstr( s >= -bigM * self.G[J,I,k+1])

        # 3: M(k-1)*A \in G(K+1)
        for i in range(1,self.m+1):
            for j in range(1,self.n+1):
                I=self.subsystem_U[i]
                J=self.subsystem_X[j]
                s=LinExpr()
                for p in range(1,self.n+1):
                    s.add(self.theta[self.K-1,i,p]*self.A[p,j]) 
                model.addConstr( s <= bigM * self.G[J,I,self.K+1])
                model.addConstr( s >= -bigM * self.G[J,I,self.K+1])                                
        
        model.addConstr( self.rho >= self.beta)
        model.addConstr( self.rho >= self.gamma)
        
        model.optimize()
        if model.Status==3:
            print "RCI is empty"
            return False
        else:
            print "\n\tGreat! RCI set exists and was computed!\n"
            for k in range(0,self.K):
                for i in range(1,self.m+1):
                    for j in range(1,self.n+1):
                        self.theta[k,i,j]=self.theta[k,i,j].X
            self.compute_D()
            for i in range(1,self.N+1):
                for j in range(1,self.N+1):
                    for k in range(1,self.K+1):
                        if i!=j:
                            self.C[k,i,j]=round(self.G[i,j,k].X)
                        else:
                            self.C[k,i,j]=1
            return True
            
            
    
    
    
    
    
                    
    # Check if a point is in RCI set         
    def is_RCI(self,x):
        w={}
        model=Model("RCIS_check")
        for k in range(0,self.K):
            for i in range(1,self.n+1):
                w[k,i]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY)
        model.update()
        # w in (1-alpha)^-1 W
        for k in range(0,self.K):
            for i in range(1,self.nW+1):
                s=LinExpr()
                for p in range(1,self.n+1):
                    s.add(self.Hw[i,p] * w[k,p])
                model.addConstr( (1-self.alpha)*s <= self.hw[i])
        # x  = xbar + D * w
        for i in range(1,self.n+1):
            s=LinExpr()
            for k in range(0,self.K):
                for p in range(1,self.n+1):
                    s.add(self.D[k,i,p] * w[k,p])
            model.addConstr(s==x[i])
        model.setParam('OutputFlag', 0) 
        model.optimize()
        if model.Status==3:
            return False
        else:
            return True
                
    def compute_control_delay(self):
        for i in range(1,self.m+1):
            self.u[i]=0
            # First term: -theta(K-1)*A*x[t-K]
            for p in range(1,self.n+1):
                for j in range(1,self.n+1):
                    self.u[i]+=-self.theta[self.K-1,i,p]*self.A[p,j]*self.xH[self.t-self.K,j]
            # Second term: [theta(k-1) -theta(k-2)A ] * x[t-k+1] from k=2 to K
            for k in range(2,self.K+1):
                for j in range(1,self.n+1):
                    self.u[i]+=self.theta[k-1,i,j]*self.xH[self.t-k+1,j]
                    for p in range(1,self.n+1):
                        self.u[i]+=-self.theta[k-2,i,p]*self.A[p,j]*self.xH[self.t-k+1,j]
            # Third term: theta(0) * x[t] from k=2 to K
            for p in range(1,self.n+1):
                self.u[i]+=self.theta[0,i,p]*self.xH[self.t,p]
            # Fourth term: - theta(k)*B*u(t-k-1) for k=0 to K-1                     
            for k in range(0,self.K):
                for p in range(1,self.n+1):
                    for j in range(1,self.m+1):
                        self.u[i]+=-self.theta[k,i,p]*self.B[p,j]*self.uH[self.t-k-1,j]
    
    
    def apply_control_delay(self,w_now):
        # Apply controls
        for i in range(1,self.m+1):
            self.uH[self.t,i]=self.u[i]
        for i in range(1,self.n+1):
            self.xH[self.t+1,i]=0
            for j in range(1,self.n+1):
                self.xH[self.t+1,i]+=self.A[i,j]*self.xH[self.t,j]
            for j in range(1,self.m+1):
                self.xH[self.t+1,i]+=self.B[i,j]*self.u[j]  
            self.xH[self.t+1,i]+=w_now[i]
        self.t+=1
        for i in range(1,self.n+1):
            self.x[i]=self.xH[self.t,i]

    def compute_control_disturbance(self):
        for i in range(1,self.m+1):
            self.u[i]=0
            for k in range(0,self.K):
                for j in range(1,self.n+1):
                    self.u[i]+=self.w[-k,j]*self.theta[k,i,j]
                        
    def evolve_disturbance(self,w_now):
        x_new={}
        for i in range(1,self.n+1):
            x_new[i]=0
            for j in range(1,self.n+1):
                x_new[i]+=self.A[i,j]*self.x[j]
            for j in range(1,self.m+1):
                x_new[i]+=self.B[i,j]*self.u[j]  
            x_new[i]+=w_now[i]
        for i in range(1,self.n+1):
            self.x[i]=x_new[i]
        w_new={}
        for k in range(1,self.K):
            for i in range(1,self.n+1):
                w_new[-k,i]=self.w[-k+1,i]
        for k in range(1,self.K):
            for i in range(1,self.n+1):
                self.w[-k,i]=w_new[-k,i] 
        for i in range(1,self.n+1):
            self.w[0,i]=w_now[i]
    

def matricize(D):
    i={}
    if isinstance( D.keys()[0], int):
        N=1
        A=np.zeros((2,1))
    else:
        N=len(D.keys()[0])
        for n in range(0,N):
            i[n]=1
        for d in D.keys():
            for n in range(0,N):
                i[n]=max(i[n],d[n])
        gen=[]        
        for n in range(0,N):
            gen.append(i[n])
    A=np.zeros(tuple(gen))
    for key in D.keys():
        s=list(key)
        s[:] = [x - 1 for x in s]
        A[tuple(s)]=D[key]
    return A