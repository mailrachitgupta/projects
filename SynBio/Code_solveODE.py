# -*- coding: utf-8 -*-
"""
Created on Mon Feb 29 22:37:24 2016

@author: Rachit
"""
import scipy as sc
import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt
import pylab
# Synthetic Biology
"""
The Following Simulations Documentations:

Basic Operations:
1. Inhibition(Symbol: '--|' ): This means in presence of an inhibitor the Product would be not formed,
   and otherwise it will be continuously formed.
   For Instance. A --| B This means that if A is present, it would stop formation of B.
2. Promotion(Symbol: '-->' ) : This means in presence of the product will be formed only when the 
   promoter is present in the vicinity.
   For Instance. A --> B This means that B will be produced only if A is present.

It is important to know that all of the reactions modelled here are simply based on principles of
chemical reaction engineering. The most beautiful aspect of these kinds of circuitry is that all of 
them have been at some point been realised in actual Biological System. Where the output is actually 
found to be similar to those of simulations. 

Below Certain major components used in Synthetic Biology have been stimulated.

"""
# Simulating an Ring Oscillator

# Basic Data


k_pdeg=0.5 #Degradation rate of Proteins (Protein degraded per unit time)
n=2 # Cooperativity Constant, this constant explains the number of units acting to reinforce
     # completion of an operation.
basal_prod=60*(5e-4) # basal_production Rate constant
k_tln= 500 #Rate of Translation
k=1 #Binding Coefficient
k_mdeg=0.2 # Reaction Constant


#Without mRNA Dynamics
def RingOscillator(v,t):
    """
        
    Published by Elowitz & Leibler. 2000. Nature
    The repressilator is a simple 3-node ring oscillator with the following structure:
    1. P1 --| P2
    2. P2 --| P3
    3. P3 --| P1
    
    In the Paper Published P1 is cI, P2 is Lac Inhibitor, P3 is TetR. The actual ciruit was realised
    in E. coli. 
    
    The below code is for simulation without considering mRNA dynamics.
    
    """
    P1,P2,P3=v
    dvdt= [basal_prod+k_tln*1/(1+(P3/k)**n)-k_pdeg*P1,
           basal_prod+k_tln*1/(1+(P1/k)**n)-k_pdeg*P2, 
           basal_prod+k_tln*1/(1+(P2/k)**n)-k_pdeg*P3,]
    return dvdt


#With mRNA Dynamics
def RingOscillatorWithmRNA(v,t):
    """
    The Logic Circuit is same as Above taken from Elowitz & Leibler. 2000. Nature
    
    The Below code takes into account mRNA dynamics. 

    Naming
    M[1,2,3] : mRNA for Protein 1,2,3 respectively    
    P[1,2,3] : Protein 1,2,3 respectively
    
    
    """    
    
    M1,P1,M2,P2,M3,P3=v #Assigning Values to the Corresponding Variables
    
    # Below Differential Equations are for the RingOscillator with mRNA kinetics. 
    dvdt= [basal_prod+k_tln*1/(1+(P3/k)**n1)-M1,beta*M1-k_pdeg*P1, 
           basal_prod+k_tln*1/(1+(P1/k)**n2)-M2,beta*M2-k_pdeg*P2, 
           basal_prod+k_tln*1/(1+(P2/k)**n3)-M3,beta*M3-k_pdeg*P3]
    
    return dvdt




#Relaxation Oscillator
def RelaxationOscillator(v,t):
    
    """
    It is yet another type of Oscillator which is used in SynBio. Involving two proteins at a time.  
    It has slightly complex behaviour. 
    This was realized by Danino et al. 2010. Nature in E. Coli in Jeff Hasty's Lab at UCSD.
    
    Model:
    1. P1 --> P2
    2. P1 --> P1
    3. P2 --| P1
    """
    
    k=1
    P1,P2=v
    dvdt=[basal_prod+100*k_tln*1/(1+(P2/k)**4)*1/(1+(10/P1)**2)-k_pdeg*P1,
          basal_prod+k_tln*1/(1+(k/P1)**4)-10*k_pdeg*P2]
    
    return dvdt

#Toggle Switch
def ToggleSwitch(v,t,toggle_action=True):
    '''
    This is a Simulation Run for Toggle Switch. In this the circuit remembers,
    the input that was given to the system. The Circuit can be simplified 
    as this
    P1 --| P2
    P2 --| P1
    a--|(P1 --| P2)
    b--|(P2 --| P1)
    '''    
    P1,P2=(v) #Protein Concentrations of P1 and P2 Species. 
    a=0 # External Inhibitor 'a'
    b=1 # External Inhibitor 'b'
    if toggle_action==True:
        if t<50:
            a=0.95
            b=0.01
        if t>50: # Switch Input After 50 time steps
            a=0.01
            b=0.95
   
    dvdt=[(basal_prod+k_tln*1/(1+(P2*(1-b)/k)**n))-k_pdeg*P1,
         (basal_prod+k_tln*1/(1+(P1*(1-a)/k)**n))-k_pdeg*P2,]
    
    return dvdt

# Boundary Conditions for The Differential Equations

TimeLength=np.linspace(0,150,200)

RingOscillator0=[1,0,0] #Initital Conditions for Ring Oscillator.
RingOscillatorWithmRNA0=[0.2,0.1,0.3,0.1,0.4,0.5] #Initital Conditions for Ring Oscillator with mRNA dynamics.
RelaxationOscillator0=[10,0] # Initial Conditions for Relaxation Oscillator.
ToggleSwitch0=[0,0] #Initial Conditions for Toggle Switch.

# Solutions of the simulations
        
SolutionTS=integrate.odeint(ToggleSwitch,ToggleSwitch0,TimeLength)
SolutionRingO=integrate.odeint(RingOscillator,RingOscillator0,TimeLength)
SolutionRelO=integrate.odeint(RelaxationOscillator,RelaxationOscillator0,TimeLength)


# Plotting of the Solutions
p1=pylab.figure(1)
pylab.plot(TimeLength,SolutionTS)
pylab.xlabel("Time")
pylab.ylabel("Concentrations")
pylab.title("Concentration versus Time Plot for Toggle Switch")

p2=pylab.figure(2)
pylab.plot(TimeLength,SolutionRingO)
pylab.xlabel("Time")
pylab.ylabel("Concentrations")
pylab.title("Concentration versus Time Plot for Ring Oscillator")

p3=pylab.figure(3)
pylab.plot(TimeLength,SolutionRelO)
pylab.xlabel("Time")
pylab.ylabel("Concentrations")
pylab.title("Concentration versus Time Plot for Relaxation Oscillator")

plt.show()

#Solving for Relaxation Oscillator
"""
for i in length:
    a.append(np.sin(i))
    pl=plt.figure(1)
    b.append(0.7)
    c.append(-0.7)

ax=plt.subplot(111)
ax.set_xlim(0,10)
ax.set_ylim(-2,2)    
plt.plot(length,a)
plt.plot(length,b)
plt.plot(length,c)
plt.show()"""
#sol = sc.optimize.odeint(rate_v, y0, t)

#def rate_p2(y,z,t):
#    return basal_prod+1/(1+(z/k)**n)-k_pdeg*y
#def rate_p3(z,x,t):
#    return basal_prod+1/(1+(x/k)**n)-k_pdeg*z


"""

#print sol
p1=plt.figure(1)
plt.plot(t, sol[:, 1])
plt.plot(t, sol[:,3])
plt.plot(t, sol[:,5])
plt.show()
"""

#print sol2
#p2=plt.figure(2)
#plt.plot(t,sol2[:,0])
#plt.plot(t,sol2[:,1])
#plt.show()

def fi(x,n):
    return (x/(30+x))**n

t=np.linspace(0,500,500)

l1=[]
l2=[]
for i in t:
    l1.append(fi(i,4))   
for i in t:
    l2.append(fi(i,2))


#p3=plt.figure(3)
#plt.plot(t,l1,'r+')
#plt.plot(t,l2)

