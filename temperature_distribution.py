#%%
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 18:41:29 2020

@author: Renan de Souza Luiz
"""
""" Método (Met):
    1 - Para o método explícito
    2 - Para o método de Euler Implícito
    3 - Para o método de Crank-Nicolson
"""

""" Item:
    a - item (a) da primeira tarefa
    b - item (b) da primeira tarefa
    c - item (c) da primeira tarefa
"""
import numpy
from matplotlib import pyplot
import math

Met=int(input("Digite o método a ser usado:"))

if Met == 1:
    
    item=str(input("Digite o item do tarefa (ex: a,b ou c):"))
    N=int(input("Digite o valor de N:"))
    lamda=float(input("Digite o valor de lambda:"))
    
    L=1
    T=1
    delx=L/N
    delt=lamda*(delx*delx)
    M=round(T/delt)
    'vetores de posição e tempo'
    x=numpy.linspace(0, 1, N + 1)
    t=numpy.linspace(0, 1, M + 1)
    To=0
    
    if item == 'a':
        
        "Ítem a)"
        
        'Condições iniciais'
        u=numpy.zeros(N+1)*To
        for i in range (0,N+1):  #criando a condição inicial
            u[i]=x[i]*x[i]*(1-x[i])*(1-x[i])
        f=numpy.zeros(N+1)
              
        'Simulação'
        num=0
        for k in range(0,M+1): 
            un=u.copy()
         
            for i in range(0,N+1):
                f[i]=10*math.cos(10*t[k])*x[i]*x[i]*(1-x[i])*(1-x[i])-(1+math.sin(10*t[k]))*(12*x[i]*x[i]-12*x[i]+2)
            
            
            if k>0:
           
                u[1:-1]=un[1:-1]+(lamda*(un[2:]-2*un[1:-1]+un[0:-2]))+delt*f[1:-1]
            
            if k == num*round(M/10):
                pyplot.plot(x,u,label=str(t[k]))
                pyplot.legend()
                pyplot.title('solução numérica - Método explícito')
                pyplot.xlabel('Pontos da barra')
                pyplot.ylabel('Temperatura')
                num=num+1
        
    
        
        
        "Solução exata"
        ue=numpy.ones(N+1)*To
        trunc=numpy.zeros(N+1)
        num=0
        plotue=pyplot.figure()
        for k in range(0,M+1): 
            
            for i in range(0,N+1):
                ue[i]=(1+math.sin(10*t[k]))*x[i]*x[i]*(1-x[i])*(1-x[i])
            for i in range(1,N):
                if k>0 and k<M:
                    trunc[i]=(((1+math.sin(10*t[k+1]))*x[i]*x[i]*(1-x[i])*(1-x[i])-ue[i])/delt)-(ue[i-1]-2*ue[i]+ue[i+1])/(delx*delx)-10*math.cos(10*t[k])*x[i]*x[i]*(1-x[i])*(1-x[i])+(1+math.sin(10*t[k]))*(12*x[i]*x[i]-12*x[i]+2)
                
            if k == num*round(M/10):
                pyplot.plot(x,ue,label=str(t[k]))
                pyplot.legend()
                pyplot.title('solução exata')
                pyplot.xlabel('Pontos da barra')
                pyplot.ylabel('Temperatura')
                num=num+1
            
        
        "Calculando o erro para T=1"
        ef=abs(u-ue)
        plotef=pyplot.figure()
        pyplot.plot(x,ef,'k')
        pyplot.title('Erro para T=1')
        pyplot.xlabel('Pontos da barra')
        pyplot.ylabel('Erro')
        print(max(ef))
        "Calculando  o erro de truncamento para T=1"
        plotef=pyplot.figure()
        pyplot.plot(x,trunc,'k')
        pyplot.title('Erro de truncamento para T=1')
        pyplot.xlabel('Pontos da barra')
        pyplot.ylabel('Erro')
    
   
    if item == 'b':
        
        "Ítem b)"
        "Condições inciais"
        ud=numpy.zeros(N+1)
        fd=numpy.zeros(N+1)
        
        for i in range(0,N+1):
            ud[i]=math.exp(-x[i])
              
        'Simulação'
        num=0
        plotud=pyplot.figure()
        for k in range(0,M+1): 
            un=ud.copy()
            ud[0]=math.exp(t[k])
            ud[-1]=math.exp(t[k]-L)*math.cos(5*t[k]*L)
            
            for i in range(0,N+1):
                fd[i]=(25*t[k]*t[k])*math.exp(t[k]-x[i])*math.cos(5*t[k]*x[i])-10*t[k]*math.exp(t[k]-x[i])*math.sin(5*t[k]*x[i])-5*x[i]*math.exp(t[k]-x[i])*math.sin(5*t[k]*x[i])
            
            if k>0:
           
                ud[1:-1]=un[1:-1]+(lamda*(un[2:]-2*un[1:-1]+un[0:-2]))+delt*fd[1:-1]
            
            if k == num*round(M/10):
                pyplot.plot(x,ud,label=str(t[k]))
                pyplot.legend()
                pyplot.title('solução numérica - Método explícito')
                pyplot.xlabel('Pontos da barra')
                pyplot.ylabel('Temperatura')
                num=num+1
            
         
        "Solução exata"
        
        ue=numpy.ones(N+1)*To
        plotuef=pyplot.figure()
        num=0
        for k in range(0,M+1): 
            
            for i in range(0,N+1):
                ue[i]=math.exp(t[k]-x[i])*math.cos(5*t[k]*x[i])
            
            if k == num*round(M/10):
                
                pyplot.plot(x,ue,label=str(t[k]))
                pyplot.legend()
                pyplot.title('solução exata')
                pyplot.xlabel('Pontos da barra')
                pyplot.ylabel('Temperatura')
                num=num+1
            
        
        "Calculando o erro para T=1"
        ef=ue-ud
        plotef=pyplot.figure()
        pyplot.plot(x,ef,'k')
        pyplot.title('Erro para T=1')
        pyplot.xlabel('Pontos da barra')
        pyplot.ylabel('Erro')
        pyplot.show(plotef)       
    
    if item == 'c':
        "Ítem c"
        h=delx
        p=0.25
        fc=numpy.zeros(N+1)
        T_num=numpy.ones(N+1)*To
        num=0
        for k in range(0,M+1): 
            
            Tnu=T_num.copy()
            for i in range(0,N+1):
                if x[i]>=p-h/2 and x[i]<=p+h/2:
                    r=10000*(1-2*t[k]*t[k])/h
                    
                    fc[i]=r
            if k>0:
                T_num[1:-1]=Tnu[1:-1]+(lamda*(Tnu[2:]-2*Tnu[1:-1]+Tnu[0:-2]))+delt*fc[1:-1]
                
            
           
          
            if t[k]>=0.1*num-delt/2 and t[k]<=0.1*num+delt/2: #k == num*round(M/10):
                pyplot.plot(x,T_num,label=str(t[k]))
                pyplot.legend()
                pyplot.title("solução numérica - Método explícito")
                pyplot.xlabel('Pontos da barra')
                pyplot.ylabel('Temperatura')
                num=num+1
                
if Met == 2:
    item=str(input("Digite o item do tarefa (ex:a,b,c):"))
    N=int(input("Digite o valor de N:"))
    L=1
    T=1
    delx=L/N
    delt=delx
    lamda=delt/(delx*delx)
    M=round(T/delt)
    

    "Criando a matriz A"
    
    ad=numpy.array([])
    al=numpy.array([])
    
    for i in range(0,N-1):
        ad=numpy.append(ad,1+2*lamda)
        al=numpy.append(al,-lamda)
        
    A=numpy.append([ad],[al],axis=0)    
    
    "Criando os vetores l e d"
    def decomposição_matrizA_LDL(ad,al,A,N,lamda):
        l=numpy.zeros(N-2)
        d=numpy.zeros(N-1)
        d[0]=ad[0]
        l[0]=al[0]/d[0]
        for i in range(N-1):
            if i<N-2:
                d[i+1]=d[0]-((l[i]*l[i])*d[i])
            if i<N-3:
                l[i+1]=-lamda/d[i+1]
        
        return [d,l]
    
    [d,l]=decomposição_matrizA_LDL(ad,al,A,N,lamda)
    
    if item == 'a':
        "criação dos vetores"
        To=0
        x=numpy.linspace(0,1,N+1)
        t=numpy.linspace(0,1,M+1)
        u=numpy.ones(N+1)*To
        for i in range (0,N+1):  #criando a condição inicial
            u[i]=x[i]*x[i]*(1-x[i])*(1-x[i])
            
        
        y=numpy.zeros(N-1)
        b=numpy.zeros(N-1)
        z=numpy.zeros(N-1)
        
        "As seguintes funções resolvem os sistemas Ly=b, Dz=y e L'u=z"
        def calcula_y(y,u,b,k,l,N,delt):
            b[0]=u[1]+delt*(10*math.cos(10*t[k])*x[1]*x[1]*(1-x[1])*(1-x[1])-(1+math.sin(10*t[k]))*(12*x[1]*x[1]-12*x[1]+2))
            y[0]=b[0]
            for i in range(0,N-2):
                b[i+1]=u[i+2]+delt*(10*math.cos(10*t[k])*x[i+2]*x[i+2]*(1-x[i+2])*(1-x[i+2])-(1+math.sin(10*t[k]))*(12*x[i+2]*x[i+2]-12*x[i+2]+2))
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
        for k in range(0,M+1):
            un=u.copy()
            yn=y.copy()
            zn=z.copy()
            
            if k>0:
                y=calcula_y(yn,un,b,k,l,N,delt)
                z=calcula_z(zn,d,y,N)
                u=calcula_u(N,z,l,un)
                u[0]=0
                
            
            if k == num*round(M/10):
                pyplot.plot(x,u,label=str(t[k]))
                pyplot.legend()
                pyplot.title('solução numérica - Euler Implícito')
                pyplot.xlabel('Pontos da barra')
                pyplot.ylabel('Temperatura')
                num=num+1
                
        "Solução exata"
        ue=numpy.ones(N+1)*To
        trunc=numpy.zeros(N+1)
        num=0
        plotue=pyplot.figure()
        for k in range(0,M+1): 
            
            for i in range(0,N+1):
                ue[i]=(1+math.sin(10*t[k]))*x[i]*x[i]*(1-x[i])*(1-x[i])
            for i in range(1,N):
                if k>0 and k<M:
                    trunc[i]=(((1+math.sin(10*t[k+1]))*x[i]*x[i]*(1-x[i])*(1-x[i])-ue[i])/delt)-(ue[i-1]-2*ue[i]+ue[i+1])/(delx*delx)-10*math.cos(10*t[k])*x[i]*x[i]*(1-x[i])*(1-x[i])+(1+math.sin(10*t[k]))*(12*x[i]*x[i]-12*x[i]+2)
                
            if k == num*round(M/10):
                pyplot.plot(x,ue,label=str(t[k]))
                pyplot.legend()
                pyplot.title('Temperatura x Pontos da barra para a solução exata')
                pyplot.xlabel('Pontos da barra')
                pyplot.ylabel('Temperatura')
                num=num+1
            
        
        "Calculando o erro para T=1"
        erro=abs(u-ue)
        plotef=pyplot.figure()
        pyplot.plot(x,erro,'k')
        pyplot.title('Erro para T=1')
        pyplot.xlabel('Pontos da barra')
        pyplot.ylabel('Erro')
        
        "Calculando  o erro de truncamento para T=1"
        plotef=pyplot.figure()
        pyplot.plot(x,trunc,'k')
        pyplot.title('Erro de truncamento para T=1')
        pyplot.xlabel('Pontos da barra')
        pyplot.ylabel('Erro')
       

    
    if item == 'b':
        "criação dos vetores"
        To=0
        x=numpy.linspace(0,1,N+1)
        t=numpy.linspace(0,1,M+1)
        u=numpy.ones(N+1)*To
        for i in range(0,N+1):
            u[i]=math.exp(-x[i])
        
        y=numpy.zeros(N-1)
        b=numpy.zeros(N-1)
        z=numpy.zeros(N-1)
        
        def calcula_y(y,u,b,k,l,N,delt,lamda):
            b[0]=u[1]+delt*((25*t[k]*t[k])*math.exp(t[k]-x[1])*math.cos(5*t[k]*x[1])-10*t[k]*math.exp(t[k]-x[1])*math.sin(5*t[k]*x[1])-5*x[1]*math.exp(t[k]-x[1])*math.sin(5*t[k]*x[1]))+lamda*math.exp(t[k])
            y[0]=b[0]
            for i in range(0,N-2):
                b[i+1]=u[i+2]+delt*((25*t[k]*t[k])*math.exp(t[k]-x[i+2])*math.cos(5*t[k]*x[i+2])-10*t[k]*math.exp(t[k]-x[i+2])*math.sin(5*t[k]*x[i+2])-5*x[i+2]*math.exp(t[k]-x[i+2])*math.sin(5*t[k]*x[i+2]))
                if i == N-3:
                    b[i+1]=b[i+1]+lamda*math.exp(t[k]-L)*math.cos(5*t[k]*L)
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
        for k in range(0,M+1):
            u[0]=math.exp(t[k])
            u[-1]=math.exp(t[k]-L)*math.cos(5*t[k]*L)
            un=u.copy()
            yn=y.copy()
            zn=z.copy()
            
            if k>0:
                y=calcula_y(yn,un,b,k,l,N,delt,lamda)
                z=calcula_z(zn,d,y,N)
                u=calcula_u(N,z,l,un)
                
                
            
            if k == num*round(M/10):
                pyplot.plot(x,u,label=str(t[k]))
                pyplot.legend()
                pyplot.title('solução numérica - Euler Implícito')
                pyplot.xlabel('Pontos da barra')
                pyplot.ylabel('Temperatura')
                num=num+1
                
       
        "Solução exata"
        
        ue=numpy.ones(N+1)*To
        plotuef=pyplot.figure()
        num=0
        for k in range(0,M+1): 
            
            for i in range(0,N+1):
                ue[i]=math.exp(t[k]-x[i])*math.cos(5*t[k]*x[i])
            
            if k == num*round(M/10):
                
                pyplot.plot(x,ue,label=str(t[k]))
                pyplot.legend()
                pyplot.title('solução exata')
                pyplot.xlabel('Pontos da barra')
                pyplot.ylabel('Temperatura')
                num=num+1
            
        
        "Calculando o erro para T=1"
        erro=abs(u-ue)
        print(max(erro))
        plotef=pyplot.figure()
        pyplot.plot(x,erro,'k')
        pyplot.title('Erro para T=1')
        pyplot.xlabel('Pontos da barra')
        pyplot.ylabel('Erro')
        
       
    if item == 'c':
        To=0
        x=numpy.linspace(0,1,N+1)
        t=numpy.linspace(0,1,M+1)
        u=numpy.ones(N+1)*To
        fc=numpy.zeros(N+1)
        y=numpy.zeros(N-1)
        b=numpy.zeros(N-1)
        z=numpy.zeros(N-1)
        
        def calcula_y(y,u,b,k,l,N,delt,lamda,fc):
            b[0]=u[1]+delt*(fc[1])
            y[0]=b[0]
            
            for i in range(0,N-2):
                b[i+1]=u[i+2]+delt*(fc[i+1])
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
        p=0.25
        num=0
        h=delx
        for k in range(0,M+1):
            un=u.copy()
            yn=y.copy()
            zn=z.copy()
            for i in range(0,N+1):
                if x[i]>=p-h/2 and x[i]<=p+h/2:
                    r=10000*(1-2*t[k]*t[k])/h
                    
                    fc[i]=r
                    
            
            if k>0:
                y=calcula_y(yn,un,b,k,l,N,delt,lamda,fc)
                z=calcula_z(zn,d,y,N)
                u=calcula_u(N,z,l,un)
                
                
            
            if t[k]>=0.1*num-delt/2 and t[k]<=0.1*num+delt/2: #k == num*round(M/10):
                pyplot.plot(x,u,label=str(t[k]))
                pyplot.legend()
                pyplot.title("solução numérica - Euler Implícito")
                pyplot.xlabel('Pontos da barra')
                pyplot.ylabel('Temperatura')
                num=num+1
        

if Met==3:
    item=str(input("Digite o item do tarefa (ex: a,b):"))
    N=int(input("Digite o valor de N:"))
    L=1
    T=1
    delx=L/N
    delt=delx
    lamda=delt/(delx*delx)
    M=round(T/delt)
   
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
    
    if item == 'a':
        "criação dos vetores"
        To=0
        x=numpy.linspace(0,1,N+1)
        t=numpy.linspace(0,1,M+1)
        u=numpy.ones(N+1)*To
        for i in range (0,N+1):  #criando a condição inicial
            u[i]=x[i]*x[i]*(1-x[i])*(1-x[i])
            
        
        y=numpy.zeros(N-1)
        b=numpy.zeros(N-1)
        z=numpy.zeros(N-1)
        
        def calcula_y(y,u,b,k,l,N,delt,lamda):
            b[0]=(1-lamda)*u[1]+(lamda/2)*(u[0]+u[2])+(delt/2)*(10*math.cos(10*t[k-1])*x[1]*x[1]*(1-x[1])*(1-x[1])-(1+math.sin(10*t[k-1]))*(12*x[1]*x[1]-12*x[1]+2)+10*math.cos(10*t[k])*x[1]*x[1]*(1-x[1])*(1-x[1])-(1+math.sin(10*t[k]))*(12*x[1]*x[1]-12*x[1]+2))
            y[0]=b[0]
            for i in range(0,N-2):
                b[i+1]=(1-lamda)*u[i+2]+(lamda/2)*(u[i+1]+u[i+3])+(delt/2)*(10*math.cos(10*t[k-1])*x[i+2]*x[i+2]*(1-x[i+2])*(1-x[i+2])-(1+math.sin(10*t[k-1]))*(12*x[i+2]*x[i+2]-12*x[i+2]+2)+10*math.cos(10*t[k])*x[i+2]*x[i+2]*(1-x[i+2])*(1-x[i+2])-(1+math.sin(10*t[k]))*(12*x[i+2]*x[i+2]-12*x[i+2]+2))
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
        for k in range(0,M+1):
            un=u.copy()
            yn=y.copy()
            zn=z.copy()
            
            if k>0:
                y=calcula_y(yn,un,b,k,l,N,delt,lamda)
                z=calcula_z(zn,d,y,N)
                u=calcula_u(N,z,l,un)
                u[0]=0
                
            
            if k == num*round(M/10):
                pyplot.plot(x,u,label=str(t[k]))
                pyplot.legend()
                pyplot.title('solução numérica - Cranck-Nicolson')
                pyplot.xlabel('Pontos da barra')
                pyplot.ylabel('Temperatura')
                num=num+1
                
                
        "Solução exata"
        
        ue=numpy.ones(N+1)*To
        trunc=numpy.zeros(N+1)
        num=0
        plotue=pyplot.figure()
        for k in range(0,M+1): 
            
            for i in range(0,N+1):
                ue[i]=(1+math.sin(10*t[k]))*x[i]*x[i]*(1-x[i])*(1-x[i])
            for i in range(1,N):
                if k>0 and k<M:
                    trunc[i]=(((1+math.sin(10*t[k+1]))*x[i]*x[i]*(1-x[i])*(1-x[i])-ue[i])/delt)-(ue[i-1]-2*ue[i]+ue[i+1])/(delx*delx)-10*math.cos(10*t[k])*x[i]*x[i]*(1-x[i])*(1-x[i])+(1+math.sin(10*t[k]))*(12*x[i]*x[i]-12*x[i]+2)
                
            if k == num*round(M/10):
                pyplot.plot(x,ue,label=str(t[k]))
                pyplot.legend()
                pyplot.title('Temperatura x Pontos da barra para a solução exata')
                pyplot.xlabel('Pontos da barra')
                pyplot.ylabel('Temperatura')
                num=num+1
            
        
        "Calculando o erro para T=1"
        erro=abs(u-ue)
        plotef=pyplot.figure()
        pyplot.plot(x,erro,'k')
        pyplot.title('Erro para T=1')
        pyplot.xlabel('Pontos da barra')
        pyplot.ylabel('Erro')
    
        "Calculando  o erro de truncamento para T=1"
        plotef=pyplot.figure()
        pyplot.plot(x,trunc,'k')
        pyplot.title('Erro de truncamento para T=1')
        pyplot.xlabel('Pontos da barra')
        pyplot.ylabel('Erro')
    
    if item == 'b':
        "criação dos vetores"
        To=0
        x=numpy.linspace(0,1,N+1)
        t=numpy.linspace(0,1,M+1)
        u=numpy.ones(N+1)*To
        for i in range(0,N+1):
            u[i]=math.exp(-x[i])
        
        y=numpy.zeros(N-1)
        b=numpy.zeros(N-1)
        z=numpy.zeros(N-1)
        
        def calcula_y(y,u,b,k,l,N,delt,lamda):
            b[0]=(1-lamda)*u[1]+(lamda/2)*(math.exp(t[k-1])+u[2])+(delt/2)*(((25*t[k-1]*t[k-1])*math.exp(t[k-1]-x[1])*math.cos(5*t[k-1]*x[1])-10*t[k-1]*math.exp(t[k-1]-x[1])*math.sin(5*t[k-1]*x[1])-5*x[1]*math.exp(t[k-1]-x[1])*math.sin(5*t[k-1]*x[1]))+((25*t[k]*t[k])*math.exp(t[k]-x[1])*math.cos(5*t[k]*x[1])-10*t[k]*math.exp(t[k]-x[1])*math.sin(5*t[k]*x[1])-5*x[1]*math.exp(t[k]-x[1])*math.sin(5*t[k]*x[1])))+(lamda/2)*(math.exp(t[k]))
            y[0]=b[0]
            for i in range(0,N-2):
                if i<N-3:
                    b[i+1]=(1-lamda)*u[i+2]+(lamda/2)*(u[i+1]+u[i+3])+(delt/2)*(((25*t[k-1]*t[k-1])*math.exp(t[k-1]-x[i+2])*math.cos(5*t[k-1]*x[i+2])-10*t[k-1]*math.exp(t[k-1]-x[i+2])*math.sin(5*t[k-1]*x[i+2])-5*x[i+2]*math.exp(t[k-1]-x[i+2])*math.sin(5*t[k-1]*x[i+2]))+((25*t[k]*t[k])*math.exp(t[k]-x[i+2])*math.cos(5*t[k]*x[i+2])-10*t[k]*math.exp(t[k]-x[i+2])*math.sin(5*t[k]*x[i+2])-5*x[i+2]*math.exp(t[k]-x[i+2])*math.sin(5*t[k]*x[i+2])))
                if i == N-3:
                    b[i+1]=(1-lamda)*u[i+2]+(lamda/2)*(u[i+1]+math.exp(t[k-1]-L)*math.cos(5*t[k-1]*L))+(delt/2)*(((25*t[k-1]*t[k-1])*math.exp(t[k-1]-x[i+2])*math.cos(5*t[k-1]*x[i+2])-10*t[k-1]*math.exp(t[k-1]-x[i+2])*math.sin(5*t[k-1]*x[i+2])-5*x[i+2]*math.exp(t[k-1]-x[i+2])*math.sin(5*t[k-1]*x[i+2]))+((25*t[k]*t[k])*math.exp(t[k]-x[i+2])*math.cos(5*t[k]*x[i+2])-10*t[k]*math.exp(t[k]-x[i+2])*math.sin(5*t[k]*x[i+2])-5*x[i+2]*math.exp(t[k]-x[i+2])*math.sin(5*t[k]*x[i+2])))+(lamda/2)*(math.exp(t[k]-L)*math.cos(5*t[k]*L))
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
        for k in range(0,M+1):
          
            un=u.copy()
            yn=y.copy()
            zn=z.copy()
            
            if k>0:
                y=calcula_y(yn,un,b,k,l,N,delt,lamda)
                z=calcula_z(zn,d,y,N)
                u=calcula_u(N,z,l,un)
                u[0]=math.exp(t[k])
                u[-1]=math.exp(t[k]-L)*math.cos(5*t[k]*L)
                
            
            if k == num*round(M/10):
                pyplot.plot(x,u,label=str(t[k]))
                pyplot.legend()
                pyplot.title('solução numérica - Cranck-Nicolson')
                pyplot.xlabel('Pontos da barra')
                pyplot.ylabel('Temperatura')
                num=num+1
                
       
        "Solução exata"
        
        ue=numpy.ones(N+1)*To
        plotuef=pyplot.figure()
        num=0
        for k in range(0,M+1): 
            
            for i in range(0,N+1):
                ue[i]=math.exp(t[k]-x[i])*math.cos(5*t[k]*x[i])
            
            if k == num*round(M/10):
                
                pyplot.plot(x,ue,label=str(t[k]))
                pyplot.legend()
                pyplot.title('solução exata')
                pyplot.xlabel('Pontos da barra')
                pyplot.ylabel('Temperatura')
                num=num+1
            
        
        "Calculando o erro para T=1"
        erro=abs(u-ue)
        plotef=pyplot.figure()
        pyplot.plot(x,erro,'k')
        pyplot.title('Erro para T=1')
        pyplot.xlabel('Pontos da barra')
        pyplot.ylabel('Erro')
    

# %%
