# -*- coding: utf-8 -*-
"""
Created on Sat Jun 18 16:48:46 2016

@author: Rachit
"""

import random
import math
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import norm
from scipy.stats import poisson

# Constants taken in the code
D=0.5 # Diffusion Constant

"""Generating Field for the movement of Bacteria and Neutrophile"""
class Field():
    def __init__(self,x,y):
        self.length=x
        self.width=y
        self.obstacle_pos=[]
        self.trace_path=[]
        self.conc_field={}
        
    def generate_obstacles(self,number=40,size=0.5):
        for i in range(number):
            self.obstacle_pos.append((random.randrange(0,self.length),random.randrange(0,self.width)))
    
    def obstacle_inPath(self,pos,cell_type):
        if cell_type=='bac':
            dia=1
        if cell_type=='neut':
            dia=5
        flag=0
        for obs_path in self.obstacle_pos:
           if calc_distance(obs_path,pos)<dia:
               flag=1
        if flag==1:
            return True
        else:
            return False

    def isInside(self,pos):
        if pos[0]<self.length and pos[1]<self.width:
            return True
        else:
            return False
    
    def store_path(self,pos):
        self.trace_path.append(pos)
    
        
    
        
    
class bacteria(object):
    def __init__(self,x,y,field):
        self.x=x
        self.y=y        
        self.field=field
        self.slength=5
        
    def get_position(self):
        return (self.x,self.y)
        
    def move(self,walls=0):
        tempx=self.x+(random.random()*2-1)
        tempy=self.y+(random.random()*2-1)
            
        if not self.field.obstacle_inPath((tempx,tempy),'bac'):
            self.x=tempx
            self.y=tempy
        else:
            self.move()
        
    def strategy_move(self,bias,pos,walls=0):
        my=-slope(pos[1],self.y)
        mx=-slope(pos[0],self.x)
       
#        rand1=math.cos(random.uniform(0,2*3.14))
#        xside=rand1*self.slength
#        yside=(1-rand1*rand1)**2*self.slength
#        rand2=random.uniform(-1,1)
#        xside_sto=rand2*self.slength
#        yside_sto=rand2*self.slength
#        tempx=self.x+(mx*bias*xside+(1-bias)*xside_sto)
#        tempy=self.y+(my*bias*yside+(1-bias)*yside_sto)
        
        rand1=random.random()
        rand2=random.random()
        rand3=random.random()
        rand4=random.random()
        
        tempx=self.x+(mx*bias*rand1+(1-bias)*(rand3*2-1))
        tempy=self.y+(my*bias*rand2+(1-bias)*(rand4*2-1))
                
        
        
        if not self.field.obstacle_inPath((tempx,tempy),'bac'):
            self.x=tempx
            self.y=tempy
        else:
            self.strategy_move(bias,pos)
    
    
class neutrophil(bacteria):      
        
    def bias_move(self,bias,slowness_factor,pos,walls=0):
        my=slope(pos[1],self.y)
        mx=slope(pos[0],self.x)
        rand1=random.random()        
        rand2=random.random()
        rand3=random.random()
        rand4=random.random()
        tempx=self.x+slowness_factor*(mx*bias*rand1+(1-bias)*(rand3*2-1))
        tempy=self.y+slowness_factor*(my*bias*rand2+(1-bias)*(rand4*2-1))
        if not self.field.obstacle_inPath((tempx,tempy),'neut'):
            self.x=tempx
            self.y=tempy
        else:
            self.bias_move(bias*0.8,slowness_factor,pos)             
    
    def grad_move(self,bias,lookforGrad=True):
        rand1=random.random()
        rand2=random.random()        
        move_dir=(rand1,rand2)
        if lookforGrad:        
            Surlist = generate_4pt((self.x,self.y))
            center_detect=conc_sum(self.field.trace_path,(self.x,self.y))
            maxGrad=0
            for i in Surlist:
                grad=conc_sum(self.field.trace_path,i)-center_detect
                if grad>maxGrad:
                    maxGrad=grad
                    move_dir=i
                    
        my=slope(move_dir[1],self.y)
        mx=slope(move_dir[0],self.x)
        rand1=random.random()        
        rand2=random.random()
        rand3=random.random()
        rand4=random.random()
        
        if maxGrad==0:
            tempx=self.x+(random.random()*2-1)
            tempy=self.y+(random.random()*2-1)
        else:
            tempx=self.x+(mx*bias*rand1+(1-bias)*(rand3*2-1))
            tempy=self.y+(my*bias*rand2+(1-bias)*(rand4*2-1))
        if not self.field.obstacle_inPath((tempx,tempy),'neut'):
            self.x=tempx
            self.y=tempy
        else:
            self.grad_move(bias*0.9,lookforGrad=True)
        
        


def conc(rCoords,t,pos=(0,0)):
    return 1/(4*3.14*t*D**2)*math.exp((-(pos[0]-rCoords[0])**2-(pos[1]-rCoords[1])**2)/D) 


""" Evaluates the concentration at any given point, provided the list containing 
the location of point sources, 'plist' is the list containing location of the
point sources and 'npos' is the location of point."""
def conc_sum(plist,npos):
    csum=0
    for t in range(len(plist[-20:])):
        csum+=conc(plist[t],len(plist[-20:])-t,npos)
    return csum

"""Generates 4 points around a given point, it is later used for getting rough idea of the 
gradient of a chemical."""
def generate_4pt(pos):
    x=5
    ptlist=[(pos[0]+x,pos[1]+x),(pos[0]-x,pos[1]+x),(pos[0]+x,pos[1]-x),\
    (pos[0]-x,pos[1]-x)]
    return ptlist
  
"""Calculates the distance between two points."""
def calc_distance(pos1,pos2):
    return math.sqrt((pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2)

def sigmoid_mod(x):
    return 1/(1+math.exp(-x))    

def bias_generator(dist):
    return 1.5-sigmoid_mod(dist)

def slope(x1,x2):
    return (x1-x2)/np.abs(x1-x2)

F1=Field(100,100)
F1.generate_obstacles(number=20)
F1_plotPoints=zip(*F1.obstacle_pos)


xstore=[]
ystore=[]
store=[]

for i in range(1):
    B1=bacteria(60,60,F1)
    N1=neutrophil(35,40,F1)
    B1_pos=B1.get_position()
    N1_pos=N1.get_position()
    B1_p=[]
    N1_p=[]
    for _ in range(1000):
        # Bacteria Movement Without any Strategy.
        #Uncomment This section to play this move.        
        bias=bias_generator(calc_distance(B1_pos,N1_pos))        
        
        
        B1.move()
        B1_pos=B1.get_position()
        
        # Bacteria Movement With Strategy #1
        # Uncomment this section to play strategic movement.
        """        
        if _%7==0:
            B1.strategy_move(bias,N1_pos)
        else:
            B1.move()
        """
        
        B1_pos=B1.get_position()
        F1.store_path(B1_pos)
        
        #Neutrophil Movement
        N1.bias_move(bias,1.2,B1_pos)
        #N1.grad_move(bias,B1_pos)
        N1_pos=N1.get_position()
        #print N1_pos,B1_pos
        
#        B1_p.append(B1_pos)
#        N1_p.append(N1_pos)        
        
        plt.clf()
        plt.xlim(0,100)
        plt.ylim(0,100)
        plt.scatter(F1_plotPoints[0],F1_plotPoints[1],c='g',s=50)
        plt.scatter(B1_pos[0],B1_pos[1],c='b',s=100)
        plt.scatter(N1_pos[0],N1_pos[1],c='r',s=8000)
        s=str(_)+".png"
        plt.savefig(s)
        if calc_distance(B1_pos,N1_pos)<=10:
            print i,"Neutrophil Caught Bacteria at step", _
            store.append(_)
            break
        if _==1000:
            store.append(1000)
        plt.pause(0.01)
        #plt.show()
        #plt.plot(bacteria1.position()[0],bacteria1.position()[1],'xb-')    


"""Plotting the Graph for Steps Distribution. """
j=sorted(store)
y=norm.pdf(j,np.mean(store),np.std(store))
plt.title("Distribution for Steps taken to catch the Bacteria")
plt.plot(j,y,'-o')
plt.xlabel("Number of Steps Taken")
plt.ylabel("Average number of cases")
plt.hist(j,normed=True)

#Some Remarks
"""
Conclusion based upon Simulations:

1. When there is random motion the spread is less, and we get narrow distribution. The bacteria
gets caught in a stipulated time with a good accuracy.
2. When there is a strategy involved there is wider spread in amount of time taken. Range is
higher for data. And we can comment with lesser confidence what time scales the bacteria
will get caught.
3. Irrespective of the field or obstacles the above two statements about spread remains 
constant. 
"""
