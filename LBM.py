#import matplotlib.pyplot as plt
import numpy
Re=100
m=100
n=100
dx=1
dt=1
e=dx/dt
l=m*dx

u=0.1
Neu=u*l/Re
Tau=3*Neu+0.5

E=[[0,0],[1,0],[0,1],[-1.0],[0,-1],[1,1],[-1,1],[-1,-1],[1,-1]]

W=[0 for i in range(9)]

W[0]=4/9
for i in range(1,5):
    W[i]=1/9

for i in range(5,9):
    W[i]=1/36

U=numpy.zeros((m,n))   
V=numpy.zeros((m,n))
Rho=numpy.ones((m,n))
Rho_U=numpy.zeros((m,n))
Rho_V=numpy.zeros((m,n))

#Making the velocities of topmost layer =0
for i in range(m):
    U[i][n-1]=0.1
    
Feq=numpy.zeros((m,n,9))
F=numpy.zeros((m,n,9))

#Finding Equilibrium distribution
for i in range(m):
    for j in range(n):
        Feq[i][j][0]=W[0]*Rho[i][j]*(1-1.5*(u**2))
        for k in range(1,9):
            Feq[i][j][k]=W[k]*Rho[i][j]*(1+3*((E[k][1])+(E[k][1]))+4.5*(((E[k][0])*U[i][j]+(E[k][1])*V[i][j])**2)-(((U[i][j])**2)+((V[i][j])**2)))

#Making initial value of F=Feq    
for i in range(m):
    for j in range(n):
        for k in range(9):
            F[i][j][k]=Feq[i][j][k]
            
#Collision step
for i in range(m):
    for j in range(n):
        for k in range(9):
            F[i][j][k]=(1-1/Tau)*F[i][j][k]+(1/Tau)*F[i][j][k]
            
#Streaming Step
for i in range(m):
    for j in range(n):
        for k in range(9):
            if E[k][0]<0:
                p=m-i
            if E[k][1]<0:
                q=n-j
            F[p][q][k]=F[p-E[k][0]][q-E[k][1]][k]

#Bounce back criteria          
for i in range(m):
    for k in [2,5,6]:
        F[i][0][k]=-F[i][0][k+2] #for 

for j in range(n):
    for k in [1,5]:
        F[0][j][k]=-F[0][j][k+2]
    F[0][j][8]=F[0][j][6]
    
    for k in [3,7]:
        F[n][j][k]=-F[n][j][k-2]
        F[n][j][6]=-F[n][j][8]
        
for i in range(m):
    for j in range(n):
        for k in range(9):
            Rho[i][j]=Rho[i][j]+F[i][j][k]
            Rho_U[i][j]=Rho_U[i][j]+F[i][j]*E[k][0]
            Rho_V[i][j]=Rho_V[i][j]+F[i][j]*E[k][1]
            
        U[i][j]=Rho_U[i][j]/Rho[i][j]
        V[i][j]=Rho_V[i][j]/Rho[i][j]
        
    U[i][n]=0
    V[i][n]=0
