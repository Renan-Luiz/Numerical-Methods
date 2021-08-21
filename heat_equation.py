# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 17:44:01 2020
@author: Renan de Souza Luiz  N°USP:9836120

Como usar o programa:
--Ao iniciar e rodar o programa deve-se escolher um dos testes:
Digite:
a para o teste 'a' do enunciado
b para o teste 'b' do enunciado 
c para o teste 'c' do enunciado
d para o teste 'd' do enunciado
--Se a escolha for o teste 'c'
Digite o valor de N
--Se a escolha for o teste 'd'
Digite o valor de N
"""
import numpy
from matplotlib import pyplot
import math
import random

Te=str(input('Digite o teste (ex:a,b,c ou d):'))

if Te=='a':

    N=128
    'numero de fontes nf'
    nf=1 
    p=[0.35]
    L=1
    T=1
    delx=L/N
    delt=delx
    lamda=delt/(delx*delx)
    M=round(T/delt)
    h=delx   
    
    "Criando a matriz A que resolve o sistema u(t,xi)"
    
    ad=numpy.array([])
    al=numpy.array([])
    
    for i in range(0,N-1):
        ad=numpy.append(ad,1+lamda)
        al=numpy.append(al,-lamda/2)
        
    A=numpy.append([ad],[al],axis=0)    
    
    def decomposição_matrizA_LDL(ad,al,A,N,lamda):
        l=numpy.zeros(N-2)
        d=numpy.zeros(N-1)
        d[0]=ad[0]
        l[0]=al[0]/d[0]
        for i in range(0,N-1):
            if i<N-2:
                d[i+1]=d[0]-((l[i]*l[i])*d[i])
            if i<N-3:    
                l[i+1]=-lamda/(2*d[i+1])
        
        return [d,l]
    
    [d,l]=decomposição_matrizA_LDL(ad,al,A,N,lamda)
    
    "criação dos vetores"
    To=0
    x=numpy.linspace(0,1,N+1)
    t=numpy.linspace(0,1,M+1)
    y=numpy.zeros(N-1)
    b=numpy.zeros(N-1)
    z=numpy.zeros(N-1)
    U=[]
    "condição inicial"
    u=numpy.ones(N+1)*To 
    
    "funções que resolvem o sistema LDLt*a=B"
    def calcula_y(y,u,b,k,l,N,delt,lamda,p,h,nf,kk):
        
        b[0]=(1-lamda)*u[1]+(lamda/2)*(u[0]+u[2])#+(delt/2)*()
        y[0]=b[0]
        
        for i in range(0,N-2):
            b[i+1]=(1-lamda)*u[i+2]+(lamda/2)*(u[i+1]+u[i+3])#+(delt/2)*()
            
            if x[i+2]>=p[kk]-h/2 and x[i+2]<=p[kk]+h/2:
                 r1=10*(1+math.cos(5*t[k-1]))
                 r2=10*(1+math.cos(5*t[k]))
                 b[i+1]=b[i+1]+(delt/2)*(r1+r2)*(1/h)
           
            y[i+1]=b[i+1]-l[i]*y[i]
        
            
        return y 
    
    def calcula_z(z,d,y,N):
        
        for i in range(0,N-1):
            z[i]=y[i]/d[i]
               
        return z
    
    def calcula_u(N,z,l,u):
        u[-2]=z[-1]
        for i in range(2,N):
            u[-i-1]=z[-i]-l[-i+1]*u[-i]
        return u
    
    "simulação"
    num=0
    for kk in range(0,nf):
        for k in range(0,M+1):
            un=u.copy()
            yn=y.copy()
            zn=z.copy()
            
            if k>0:
                y=calcula_y(yn,un,b,k,l,N,delt,lamda,p,h,nf,kk)
                z=calcula_z(zn,d,y,N)
                u=calcula_u(N,z,l,un)
                
            if k == 128:
                pyplot.plot(x,u,label='u'+str(kk+1)+'(T,xi)')
                pyplot.legend()
                pyplot.title('solução numérica - Cranck-Nicolson')
                pyplot.xlabel('Pontos da barra')
                pyplot.ylabel('Temperatura')
        
        U.append(u)
              
            
        
            
            
    
    "Os valores já conhecidos de u(T,xi)"
    uT=numpy.ones(N+1)*To
    for i in range(0,N+1):
        uT[i]=7*U[0][i]
    pyplot.plot(x,uT,label='uT')
    pyplot.legend()
    pyplot.title('Solução T=1')
    pyplot.xlabel('Pontos da barra')
    pyplot.ylabel('Temperatura')
       
    
    'Criando as matrizes A e B do sistema Au=B'
    A=[]
    B=[]
    for k in range(0,nf):
        ut_prod=0
        u_p=[]
        for j in range (0,nf):
            u_prod=0
          
            for i in range(1,N): 
                u_prod=u_prod+U[j][i]*U[k][i]
                if j==0:
                    ut_prod=ut_prod+uT[i]*U[k][i]
            u_p.append(u_prod)
        
        A.append(u_p)   
        B.append(ut_prod)   
           
        
    
    l=numpy.identity(nf)
    d=numpy.zeros(nf)
    v=numpy.zeros(nf)
    
    'Função que faz a decomposição A em LDLt'
    def decomposição_matrizA_LDLt(l,d,v,A,nf):
    
    
        d[0]=A[0][0]
        for i in range(0,nf):
            for j in range(0,i):
                v[j]=l[i][j]*d[j]
                #if j>i:
                    #d[i]=A[i][i]
                
                soma=0
                if j<i:
                   
                    for k in range(0,i):
                        soma=soma+l[i][k]*v[k]
                        
                d[i]=A[i][i]-soma
                 
            
            for j in range(i+1,nf):
                soma=0
                for k in range(0,i):
                    soma=soma+(l[j][k]*v[k])
                
                
                l[j][i]=(A[j][i]-soma)/d[i]
                       
                    
                    
        return [d,l]
    
    [d,l]=decomposição_matrizA_LDLt(l,d,v,A,nf)
    
    "Funções que resolvem o sistema Lz=B, Dy=x e Lta=y, respectivamente"
    def calcula_y2(l,B,nf):
        y=numpy.zeros(nf)
        y[0]=B[0]
        for i in range(1,nf):
            soma=0
            for j in range(1,i+1):
                soma=soma+l[i][j-1]*y[j-1]
            
            y[i]=B[i]-soma
        
            
        return y 
    
    
    def calcula_z2(d,y,nf):
        z=numpy.zeros(nf)
        for i in range(0,nf):
            z[i]=y[i]/d[i]
               
        return z
    
    def calcula_a2(l,z,nf):
        a=numpy.zeros(nf)
        for i in range(1,nf+1):
            soma=0
            for j in range(1,i):
                soma=soma+l[nf-j][nf-i]*a[nf-j]
            
            a[nf-i]=z[nf-i]-soma
            
            
            
        return a
    
    y2=calcula_y2(l,B,nf)
    z2=calcula_z2(d,y2,nf)
    a=calcula_a2(l,z2,nf)
    print("Os valores das intensidades a's são:"+str(a))
    
    
    
    
    soma=0
    for i in range(1,N):
        au=0
        for kk in range(0,nf):
            au=au+a[kk]*U[kk][i]
        soma=soma+(uT[i]-au)*(uT[i]-au)
    E=math.sqrt(delx*soma)
    
    print("Os valor do erro é:"+str(E))

if Te=='b':
    N=128
    nf=4
    p=[0.15,0.3,0.7,0.8]
    L=1
    T=1
    delx=L/N
    delt=delx
    lamda=delt/(delx*delx)
    M=round(T/delt)
    h=delx   
    "Criando a matriz A"
    
    ad=numpy.array([])
    al=numpy.array([])
    
    for i in range(0,N-1):
        ad=numpy.append(ad,1+lamda)
        al=numpy.append(al,-lamda/2)
        
    A=numpy.append([ad],[al],axis=0)    
    
    def decomposição_matrizA_LDL(ad,al,A,N,lamda):
        l=numpy.zeros(N-2)
        d=numpy.zeros(N-1)
        d[0]=ad[0]
        l[0]=al[0]/d[0]
        for i in range(0,N-1):
            if i<N-2:
                d[i+1]=d[0]-((l[i]*l[i])*d[i])
            if i<N-3:    
                l[i+1]=-lamda/(2*d[i+1])
        
        return [d,l]
    
    [d,l]=decomposição_matrizA_LDL(ad,al,A,N,lamda)
    
    "criação dos vetores"
    To=0
    x=numpy.linspace(0,1,N+1)
    t=numpy.linspace(0,1,M+1)
    y=numpy.zeros(N-1)
    b=numpy.zeros(N-1)
    z=numpy.zeros(N-1)
    U=[]
    "condição inicial"
    u=numpy.ones(N+1)*To 
    
    "funções que resolvem o sistema LDLt*a=B"
    def calcula_y(y,u,b,k,l,N,delt,lamda,p,h,nf,kk):
        
        b[0]=(1-lamda)*u[1]+(lamda/2)*(u[0]+u[2])#+(delt/2)*()
        y[0]=b[0]
        
        for i in range(0,N-2):
            b[i+1]=(1-lamda)*u[i+2]+(lamda/2)*(u[i+1]+u[i+3])#+(delt/2)*()
            
            if x[i+2]>=p[kk]-h/2 and x[i+2]<=p[kk]+h/2:
                 r1=10*(1+math.cos(5*t[k-1]))
                 r2=10*(1+math.cos(5*t[k]))
                 b[i+1]=b[i+1]+(delt/2)*(r1+r2)*(1/h)
           
            y[i+1]=b[i+1]-l[i]*y[i]
        
            
        return y 
    
    def calcula_z(z,d,y,N):
        
        for i in range(0,N-1):
            z[i]=y[i]/d[i]
               
        return z
    
    def calcula_u(N,z,l,u):
        u[-2]=z[-1]
        for i in range(2,N):
            u[-i-1]=z[-i]-l[-i+1]*u[-i]
        return u
    
    "simulação"
    num=0
    for kk in range(0,nf):
        for k in range(0,M+1):
            un=u.copy()
            yn=y.copy()
            zn=z.copy()
            
            if k>0:
                y=calcula_y(yn,un,b,k,l,N,delt,lamda,p,h,nf,kk)
                z=calcula_z(zn,d,y,N)
                u=calcula_u(N,z,l,un)
                
            if k == 128:
                pyplot.plot(x,u,label='u'+str(kk+1)+'(T,xi)')
                pyplot.legend()
                pyplot.title('solução numérica - Cranck-Nicolson')
                pyplot.xlabel('Pontos da barra')
                pyplot.ylabel('Temperatura')
        
        U.append(u)
              
            
        
            
            
    
    "Os valores já conhecidos de u(T,xi)"
    uT=numpy.ones(N+1)*To
    for i in range(0,N+1):
    #    uT[i]=7*U[0][i]
        uT[i]=2.3*U[0][i]+3.7*U[1][i]+0.3*U[2][i]+4.2*U[3][i]
    pyplot.plot(x,uT,label='uT')
    pyplot.legend()
    pyplot.title('Solução em T=1')
    pyplot.xlabel('Pontos da barra')
    pyplot.ylabel('Temperatura')
       
    
    
    A=[]
    B=[]
    for k in range(0,nf):
        ut_prod=0
        u_p=[]
        for j in range (0,nf):
            u_prod=0
          
            for i in range(1,N): 
                u_prod=u_prod+U[j][i]*U[k][i]
                if j==0:
                    ut_prod=ut_prod+uT[i]*U[k][i]
            u_p.append(u_prod)
        
        A.append(u_p)   
        B.append(ut_prod)   
           
        
    
    l=numpy.identity(nf)
    d=numpy.zeros(nf)
    v=numpy.zeros(nf)
    def decomposição_matrizA_LDLt(l,d,v,A,nf):
    
    
        d[0]=A[0][0]
        for i in range(0,nf):
            for j in range(0,i):
                v[j]=l[i][j]*d[j]
                #if j>i:
                    #d[i]=A[i][i]
                
                soma=0
                if j<i:
                   
                    for k in range(0,i):
                        soma=soma+l[i][k]*v[k]
                        
                d[i]=A[i][i]-soma
                 
            
            for j in range(i+1,nf):
                soma=0
                for k in range(0,i):
                    soma=soma+(l[j][k]*v[k])
                
                
                l[j][i]=(A[j][i]-soma)/d[i]
                       
                    
                    
        return [d,l]
    
    [d,l]=decomposição_matrizA_LDLt(l,d,v,A,nf)
    
    def calcula_y2(l,B,nf):
        y=numpy.zeros(nf)
        y[0]=B[0]
        for i in range(1,nf):
            soma=0
            for j in range(1,i+1):
                soma=soma+l[i][j-1]*y[j-1]
            
            y[i]=B[i]-soma
        
            
        return y 
    
    def calcula_z2(d,y,nf):
        z=numpy.zeros(nf)
        for i in range(0,nf):
            z[i]=y[i]/d[i]
               
        return z
    
    def calcula_a2(l,z,nf):
        a=numpy.zeros(nf)
        for i in range(1,nf+1):
            soma=0
            for j in range(1,i):
                soma=soma+l[nf-j][nf-i]*a[nf-j]
            
            a[nf-i]=z[nf-i]-soma
            
            
            
        return a
    
    y2=calcula_y2(l,B,nf)
    z2=calcula_z2(d,y2,nf)
    a=calcula_a2(l,z2,nf)
    print("Os valores das intensidades a's são:"+str(a))
    
    
    
    
    soma=0
    for i in range(1,N):
        au=0
        for kk in range(0,nf):
            au=au+a[kk]*U[kk][i]
        soma=soma+(uT[i]-au)*(uT[i]-au)
    E=math.sqrt(delx*soma)
    
    print("Os valor do erro é:"+str(E))
    
if Te=='c':
    N=int(input('Digite o valor de N:'))
    nf=10
    p=[0.14999999999999999,0.20000000000000001,0.29999999999999999,0.34999999999999998,0.50000000000000000,0.59999999999999998,0.69999999999999996,0.72999999999999998,0.84999999999999998,0.90000000000000002]
    L=1
    T=1
    delx=L/N
    delt=delx
    lamda=delt/(delx*delx)
    M=round(T/delt)
    h=delx   
    "Criando a matriz A"
    
    ad=numpy.array([])
    al=numpy.array([])
    
    for i in range(0,N-1):
        ad=numpy.append(ad,1+lamda)
        al=numpy.append(al,-lamda/2)
        
    A=numpy.append([ad],[al],axis=0)    
    
    def decomposição_matrizA_LDL(ad,al,A,N,lamda):
        l=numpy.zeros(N-2)
        d=numpy.zeros(N-1)
        d[0]=ad[0]
        l[0]=al[0]/d[0]
        for i in range(0,N-1):
            if i<N-2:
                d[i+1]=d[0]-((l[i]*l[i])*d[i])
            if i<N-3:    
                l[i+1]=-lamda/(2*d[i+1])
        
        return [d,l]
    
    [d,l]=decomposição_matrizA_LDL(ad,al,A,N,lamda)
    
    "criação dos vetores"
    To=0
    x=numpy.linspace(0,1,N+1)
    t=numpy.linspace(0,1,M+1)
    y=numpy.zeros(N-1)
    b=numpy.zeros(N-1)
    z=numpy.zeros(N-1)
    U=[]
    "condição inicial"
    u=numpy.ones(N+1)*To 
        
    def calcula_y1(y,u,b,k,l,N,delt,lamda,p,h,nf,kk):
        
        b[0]=(1-lamda)*u[1]+(lamda/2)*(u[0]+u[2])#+(delt/2)*()
        y[0]=b[0]
        
        for i in range(0,N-2):
            b[i+1]=(1-lamda)*u[i+2]+(lamda/2)*(u[i+1]+u[i+3])#+(delt/2)*()
            
            if x[i+2]>=p[kk]-h/2 and x[i+2]<=p[kk]+h/2:
                 r1=10*(1+math.cos(5*t[k-1]))
                 r2=10*(1+math.cos(5*t[k]))
                 b[i+1]=b[i+1]+(delt/2)*(r1+r2)*(1/h)
           
            y[i+1]=b[i+1]-l[i]*y[i]
        
            
        return y 
    
    def calcula_z1(z,d,y,N):
        
        for i in range(0,N-1):
            z[i]=y[i]/d[i]
               
        return z
    
    def calcula_u1(N,z,l,u):
        u[-2]=z[-1]
        for i in range(2,N):
            u[-i-1]=z[-i]-l[-i+1]*u[-i]
        return u
    
    "simulação"
    num=0
    for kk in range(0,nf):
        for k in range(0,M+1):
            un=u.copy()
            yn=y.copy()
            zn=z.copy()
            
            if k>0:
                y=calcula_y1(yn,un,b,k,l,N,delt,lamda,p,h,nf,kk)
                z=calcula_z1(zn,d,y,N)
                u=calcula_u1(N,z,l,un)
                
            if k == 128:
                pyplot.plot(x,u,label='u'+str(kk+1)+'(T,xi)')
                pyplot.legend()
                pyplot.title('solução numérica - Cranck-Nicolson')
                pyplot.xlabel('Pontos da barra')
                pyplot.ylabel('Temperatura')
        
        U.append(u)
              
            
        
            
            
    uTarq=open('valores de uT.txt','r')
    
    num1=0
    num2=2048/N
    uT=[]
    for i in range(0,2049):
        n=float(uTarq.readline())
        if i==num1*num2:
            uT.append(n)
            num1=num1+1
        
      
    
    uTarq.close()
    
    "Os valores já conhecidos de u(T,xi)"
    pyplot.plot(x,uT,label='uT')
    pyplot.legend()
    pyplot.title('Solução em T=1')
    pyplot.xlabel('Pontos da barra')
    pyplot.ylabel('Temperatura')
       
    
    
    A=[]
    B=[]
    for k in range(0,nf):
        ut_prod=0
        u_p=[]
        for j in range (0,nf):
            u_prod=0
          
            for i in range(1,N): 
                u_prod=u_prod+U[j][i]*U[k][i]
                if j==0:
                    ut_prod=ut_prod+uT[i]*U[k][i]
            u_p.append(u_prod)
        
        A.append(u_p)   
        B.append(ut_prod)
    
        
    
    #a=numpy.linalg.solve(A,B)
    #print("Os valores das intensidades a's são:"+str(a))
    
    l=numpy.identity(nf)
    d=numpy.zeros(nf)
    v=numpy.zeros(nf)
    def decomposição_matrizA_LDLt(l,d,v,A,nf):
    
    
        d[0]=A[0][0]
        for i in range(0,nf):
            for j in range(0,i):
                v[j]=l[i][j]*d[j]
                #if j>i:
                    #d[i]=A[i][i]
                
                soma=0
                if j<i:
                   
                    for k in range(0,i):
                        soma=soma+l[i][k]*v[k]
                        
                d[i]=A[i][i]-soma
            
            
            for j in range(i+1,nf):
                soma=0
                for k in range(0,i):
                    soma=soma+(l[j][k]*v[k])
                
                
                l[j][i]=(A[j][i]-soma)/d[i]
                       
                    
                    
        return [d,l]
    
    [d,l]=decomposição_matrizA_LDLt(l,d,v,A,nf)
    
    def calcula_y2(l,B,nf):
        y=numpy.zeros(nf)
        y[0]=B[0]
        for i in range(1,nf):
            soma=0
            for j in range(1,i+1):
                soma=soma+l[i][j-1]*y[j-1]
            
            y[i]=B[i]-soma
        
            
        return y 
    
    def calcula_z2(d,y,nf):
        z=numpy.zeros(nf)
        for i in range(0,nf):
            z[i]=y[i]/d[i]
               
        return z
    
    def calcula_a2(l,z,nf):
        a=numpy.zeros(nf)
        for i in range(1,nf+1):
            soma=0
            for j in range(1,i):
                soma=soma+l[nf-j][nf-i]*a[nf-j]
            
            a[nf-i]=z[nf-i]-soma
            
            
            
        return a
    
    y=calcula_y2(l,B,nf)
    z=calcula_z2(d,y,nf)
    a=calcula_a2(l,z,nf)
    print("Os valores das intensidades a's são:"+str(a))
    
    
    soma=0
    for i in range(1,N):
        au=0
        for kk in range(0,nf):
            au=au+a[kk]*U[kk][i]
        soma=soma+(uT[i]-au)*(uT[i]-au)
    E=math.sqrt(delx*soma)        
    
    print("Os valor do erro é:"+str(E))

if Te=='d':
    N=int(input('Digite o valor de N:'))
    nf=10
    p=[0.14999999999999999,0.20000000000000001,0.29999999999999999,0.34999999999999998,0.50000000000000000,0.59999999999999998,0.69999999999999996,0.72999999999999998,0.84999999999999998,0.90000000000000002]
    L=1
    T=1
    delx=L/N
    delt=delx
    lamda=delt/(delx*delx)
    M=round(T/delt)
    h=delx   
    e=0.01
    "Criando a matriz A"
    
    ad=numpy.array([])
    al=numpy.array([])
    
    for i in range(0,N-1):
        ad=numpy.append(ad,1+lamda)
        al=numpy.append(al,-lamda/2)
        
    A=numpy.append([ad],[al],axis=0)    
    
    def decomposição_matrizA_LDL(ad,al,A,N,lamda):
        l=numpy.zeros(N-2)
        d=numpy.zeros(N-1)
        d[0]=ad[0]
        l[0]=al[0]/d[0]
        for i in range(0,N-1):
            if i<N-2:
                d[i+1]=d[0]-((l[i]*l[i])*d[i])
            if i<N-3:    
                l[i+1]=-lamda/(2*d[i+1])
        
        return [d,l]
    
    [d,l]=decomposição_matrizA_LDL(ad,al,A,N,lamda)
    
    "criação dos vetores"
    To=0
    x=numpy.linspace(0,1,N+1)
    t=numpy.linspace(0,1,M+1)
    y=numpy.zeros(N-1)
    b=numpy.zeros(N-1)
    z=numpy.zeros(N-1)
    U=[]
    "condição inicial"
    u=numpy.ones(N+1)*To 
        
    def calcula_y1(y,u,b,k,l,N,delt,lamda,p,h,nf,kk):
        
        b[0]=(1-lamda)*u[1]+(lamda/2)*(u[0]+u[2])#+(delt/2)*()
        y[0]=b[0]
        
        for i in range(0,N-2):
            b[i+1]=(1-lamda)*u[i+2]+(lamda/2)*(u[i+1]+u[i+3])#+(delt/2)*()
            
            if x[i+2]>=p[kk]-h/2 and x[i+2]<=p[kk]+h/2:
                 r1=10*(1+math.cos(5*t[k-1]))
                 r2=10*(1+math.cos(5*t[k]))
                 b[i+1]=b[i+1]+(delt/2)*(r1+r2)*(1/h)
           
            y[i+1]=b[i+1]-l[i]*y[i]
        
            
        return y 
    
    def calcula_z1(z,d,y,N):
        
        for i in range(0,N-1):
            z[i]=y[i]/d[i]
               
        return z
    
    def calcula_u1(N,z,l,u):
        u[-2]=z[-1]
        for i in range(2,N):
            u[-i-1]=z[-i]-l[-i+1]*u[-i]
        return u
    
    "simulação"
    num=0
    for kk in range(0,nf):
        for k in range(0,M+1):
            un=u.copy()
            yn=y.copy()
            zn=z.copy()
            
            if k>0:
                y=calcula_y1(yn,un,b,k,l,N,delt,lamda,p,h,nf,kk)
                z=calcula_z1(zn,d,y,N)
                u=calcula_u1(N,z,l,un)
                
            if k == 128:
                pyplot.plot(x,u,label='u'+str(kk+1)+'(T,xi)')
                pyplot.legend()
                pyplot.title('solução numérica - Cranck-Nicolson')
                pyplot.xlabel('Pontos da barra')
                pyplot.ylabel('Temperatura')
        
        U.append(u)
              
            
        
            
            
    uTarq=open('valores de uT.txt','r')
    
    num1=0
    num2=2048/N
    uT=[]
    for i in range(0,2049):
        r=(random.random()-0.5)*2
        n_org=float(uTarq.readline())
        n=n_org*(1+r*e)
        if i==num1*num2:
            uT.append(n)
            num1=num1+1
        
      
    
    uTarq.close()
    
    "Os valores já conhecidos de u(T,xi)"
    pyplot.plot(x,uT,label='uT')
    pyplot.legend()
    pyplot.title('Solução exata em T=1')
    pyplot.xlabel('Pontos da barra')
    pyplot.ylabel('Temperatura')
       
    
    
    A=[]
    B=[]
    for k in range(0,nf):
        ut_prod=0
        u_p=[]
        for j in range (0,nf):
            u_prod=0
          
            for i in range(1,N): 
                u_prod=u_prod+U[j][i]*U[k][i]
                if j==0:
                    ut_prod=ut_prod+uT[i]*U[k][i]
            u_p.append(u_prod)
        
        A.append(u_p)   
        B.append(ut_prod)
    
        
    
    
    
    l=numpy.identity(nf)
    d=numpy.zeros(nf)
    v=numpy.zeros(nf)
    def decomposição_matrizA_LDLt(l,d,v,A,nf):
    
    
        d[0]=A[0][0]
        for i in range(0,nf):
            for j in range(0,i):
                v[j]=l[i][j]*d[j]
                #if j>i:
                    #d[i]=A[i][i]
                
                soma=0
                if j<i:
                   
                    for k in range(0,i):
                        soma=soma+l[i][k]*v[k]
                        
                d[i]=A[i][i]-soma
            
            
            for j in range(i+1,nf):
                soma=0
                for k in range(0,i):
                    soma=soma+(l[j][k]*v[k])
                
                
                l[j][i]=(A[j][i]-soma)/d[i]
                       
                    
                    
        return [d,l]
    
    [d,l]=decomposição_matrizA_LDLt(l,d,v,A,nf)
    
    def calcula_y2(l,B,nf):
        y=numpy.zeros(nf)
        y[0]=B[0]
        for i in range(1,nf):
            soma=0
            for j in range(1,i+1):
                soma=soma+l[i][j-1]*y[j-1]
            
            y[i]=B[i]-soma
        
            
        return y 
    
    def calcula_z2(d,y,nf):
        z=numpy.zeros(nf)
        for i in range(0,nf):
            z[i]=y[i]/d[i]
               
        return z
    
    def calcula_a2(l,z,nf):
        a=numpy.zeros(nf)
        for i in range(1,nf+1):
            soma=0
            for j in range(1,i):
                soma=soma+l[nf-j][nf-i]*a[nf-j]
            
            a[nf-i]=z[nf-i]-soma
            
            
            
        return a
    
    y=calcula_y2(l,B,nf)
    z=calcula_z2(d,y,nf)
    a=calcula_a2(l,z,nf)
    print("Os valores das intensidades a's são:"+str(a))
    
    
    soma=0
    for i in range(1,N):
        au=0
        for kk in range(0,nf):
            au=au+a[kk]*U[kk][i]
        soma=soma+(uT[i]-au)*(uT[i]-au)
    E=math.sqrt(delx*soma)        
    
    print("Os valor do erro é:"+str(E))

    