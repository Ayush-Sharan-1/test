# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 01:10:15 2020

@author: ayush
"""
#asdfghjk
import matplotlib.pyplot as plt
import numpy
#import matplotlib.gridspec as gridspec
Re=100
#m=100
#n=100
m=10
n=10
dx=1
dt=1
l=m*dx

u=0.1
Neu=u*l/Re
Tau=3*Neu+0.5

E=[[0,0],[1,0],[0,1],[-1,0],[0,-1],[1,1],[-1,1],[-1,-1],[1,-1]]

W=[0 for i in range(9)]

W[0]=4/9
for i in range(1,5):
    W[i]=1/9

for i in range(5,9):
    W[i]=1/36

U=numpy.zeros((m,n))   
V=numpy.zeros((m,n))
U2=numpy.zeros((m,n))   
V2=numpy.zeros((m,n))
Rho=numpy.ones((m,n))
Rho_U=numpy.zeros((m,n))
Rho_V=numpy.zeros((m,n))

#Making the x direction velocity of topmost layer =0.1
for i in range(m):
    U[i][n-1]=0.1
    
Feq=numpy.zeros((m,n,9))
F=numpy.zeros((m,n,9))

#Finding Equilibrium distribution
for i in range(m):
    for j in range(n):
       for k in range(9):
           Feq[i][j][k]=W[k]*Rho[i][j]*(1+3*((E[k][0]*U[i][j])+(E[k][1])*V[i][j])+4.5*(((E[k][0])*U[i][j]+(E[k][1])*V[i][j])**2)-1.5*(((U[i][j])**2)+((V[i][j])**2)))

z=0      
#Making initial value of F=Feq    
F=Feq.copy()

not_convergent=True
while not_convergent:          
   
    F=(1-1/Tau)*F+(1/Tau)*Feq
                
    #Streaming Step
    p=0
    q=0
    for i in range(m-1):
        for j in range(n-1):
            for k in range(9):
                if E[k][0]>0:
                    #print("1")
                    p=m-i-1
                else:
                    #print("2")
                    p=i
                if E[k][1]>0:
                    #print("3")
                    q=n-j-1
                else:
                    #print("4")
                #print(".....")
                    q=j
                F[p][q][k]=F[p-E[k][0]][q-E[k][1]][k]
    for i in range(m):
        for k in range(9):
            F[i][n-1][k]=Feq[i][n-1][k]
    
    #Bounce back criteria          
    for i in range(m):
        for k in [2,5,6]:
            F[i][0][k]=-F[i][0][k+2] 
    
    for j in range(n):
        for k in [1,5]:
            F[0][j][k]=-F[0][j][k+2]
        F[0][j][8]=F[0][j][6]
        
        for k in [3,7]:
            F[n-1][j][k]=-F[n-1][j][k-2]
            F[n-1][j][6]=-F[n-1][j][8]
            
    for i in range(m):
        for j in range(n):
            for k in range(9):
                Rho[i][j]=Rho[i][j]+F[i][j][k]
                Rho_U[i][j]+=F[i][j][k]*E[k][0]
                Rho_V[i][j]+=F[i][j][k]*E[k][1]
        
            U2[i][j]=Rho_U[i][j]/Rho[i][j]
            V2[i][j]=Rho_V[i][j]/Rho[i][j]
            
        U2[i][n-1]=0.1
        V2[i][n-1]=0
        
    #checking for convergence
    not_convergent=False
    for i in range(m-1):
        for j in range(n-1):
            if abs(U[i][j]-U2[i][j])>1E-6: #Changed to 1E-6 from 1E-8 for debugging to get the result faster
                not_convergent=True
                #print(abs(U[i][j]-U2[i][j]))
                break
    
    U=U2.copy()
    V=V2.copy()    
    
                
    Feq[i][j][k]=W[k]*Rho[i][j]*(1+3*((E[k][1])+(E[1][1]))+4.5*(((E[k][0])*U[i][j]+(E[k][1])*V[i][j])**2)-1.5*(((U[i][j])**2)+((V[i][j])**2)))
    print(z)
    z+=1
  

#x=numpy.linspace(0,99,100)
#y=numpy.linspace(0,99,100)
x=numpy.linspace(0,9,10)
y=numpy.linspace(0,9,10)
X,Y =numpy.meshgrid(x,y)

u =U.transpose()
v =V.transpose()

vels=(U**2+V**2)**0.5

lw=4*vels/numpy.max(vels)

plt.figure()
stream=plt.streamplot(X,Y,u,v,arrowsize=1,arrowstyle='->',color=vels,density=1,cmap=plt.cm.jet)

plt.title("Lid driven cavity flow")
plt.xlabel("x")
plt.ylabel("y")
plt.colorbar(stream.lines)

plt.show() 