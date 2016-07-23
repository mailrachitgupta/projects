# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 12:36:11 2016
@author: Rachit
"""

"""
Problem Statement: To Create a 4 level equivalent of Game of Life, using 4 elements:
1. Bare Earth referred as 0
2. Grass referred as 1
3. Prey referred as 2
4. Predators referred as 3
with the rules analogous to Game of Life. 

Simplifications: The border elements haven't been modified taken for simplicity.
"""
import numpy as np
import collections
import matplotlib.pyplot as plt
from matplotlib import colors

np.random.seed(1) #Setting Random Seed=1, for reproducibility(optional)

"""Start Of Class GoL_4level()"""

class GoL_4level():
    
    Z=np.random.randint(4,size=(80,80));
    
    # Setting the borders as zero    
    
    Z[0:20,:]=0
    Z[:,0:20]=0
    Z[60:-1,:]=0
    Z[:,60:-1]=0
    n=80
    
    """ Function to Search for Neighbours and return the neighbour."""
    
    def search_neigh(self,Z,num,xindex,yindex):
        search_m=Z[xindex-1:xindex+2,yindex-1:yindex+2]
        ones=0
        bin_matrix=search_m==num
        for x in range(3):
            for y in range(3):
                if bin_matrix[x,y]==True:
                    ones+=1
        if bin_matrix[1][1]==True:
            return ones-1
        else:
            return ones

    """Function to generate Neighbourhood of each element it will create a 3D array for storage."""
    
    def neighb_count(self):
        N=np.zeros((4,self.n,self.n),dtype=int)
        for x in range(1,self.n-1):
            for y in range(1,self.n-1):
                for num in range(4):                
                    N[num,x,y]=self.search_neigh(self.Z,num,x,y)
        return N
    """#!"""

    """Individual Functions for each Level"""    
    
    def when_BareEarth(self,N,x,y):
        if N[1,x,y]>=2 and N[2,x,y]>=1:
            self.Z[x,y]=2
        elif N[1,x,y]>0:
            self.Z[x,y]=1
    
    def when_Grass(self,N,x,y):
        if N[2,x,y]>2 and N[3,x,y]>=1:
            self.Z[x,y]=3
        elif N[2,x,y]>1:
            self.Z[x,y]=0
        
    def when_Prey(self,N,x,y):
        if N[1,x,y]<2:
            self.Z[x,y]=0
        if N[2,x,y]>2 and N[3,x,y]>1:
            self.Z[x,y]=3
            
    def when_Predator(self,N,x,y):
        if N[2,x,y]<2:
            self.Z[x,y]=0
        
    
    """#!"""
        
    """Function to Analyse the neighbourhood and make changes, The engine """
    
    def Analysis(self,N):
                
        for x in range(0,self.n):
            for y in range(0,self.n):
                if self.Z[x,y]==0:
                    self.when_BareEarth(N,x,y)
                elif self.Z[x,y]==1:
                    self.when_Grass(N,x,y)
                elif self.Z[x,y]==2:
                    self.when_Prey(N,x,y)
                elif self.Z[x,y]==3:
                    self.when_Predator(N,x,y)
        return self.Z 
        """ #! """
    """#! End of Class GoL_4level()"""

#Declaration of Object 

a2=GoL_4level()

# Make a Color Map of Fixed Colors

cmap = colors.ListedColormap(['sandybrown', 'green', 'yellow','red'])
bounds=[0,1,2,3]
norm = colors.BoundaryNorm(bounds, cmap.N)

"""Animations for the forest"""

for _ in range(50):
    kool= a2.neighb_count()
    img=a2.Analysis(kool)
    plt.pause(0.05)
    plt.title("Story of a Forest")
    imgplot= plt.imshow(img,vmin=0,vmax=3,cmap=cmap,interpolation='nearest')
plt.show()

"""-------------End of Story--------------"""
