"""
Created on Thu Oct 12 15:13:35 2023

Get points of line

Create baseline using 2 methods 
    1. Linear regression 
    2. Gaussian kernel 
"""
import os 
import numpy as np
import Module.Functions as fb
import shapely
from scipy import signal
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt
import pandas as pd

#===========================FUNCTIONS=================================


def findDistancesGauss(x,y,xscipy, yscipy, dx, dy):
    #variables 
    distanceE = []
    saveIndex = []
    
    realPair = list(zip(x,y))
    polyGon = shapely.geometry.LineString(realPair)
    
    for j in range(1,len(dx)):
        xs, ys = fb.createNormalLine(xscipy[j], yscipy[j], dx[j], dy[j])
       
        
        stack = np.stack((xs,ys), axis=-1)
        line = shapely.geometry.LineString(stack)
        
        #TODO remove this from main CODE
        if(polyGon.intersects(line)):
            #intersection geometry
            interPoints = polyGon.intersection(line)
            
            #intersection point
            mx, my = fb.proccessIntersectionPoint(interPoints, xscipy[j], yscipy[j])
            
            euD = fb.euclidDist(xscipy[j], yscipy[j], mx, my)
            distanceE.append(euD)
            saveIndex.append(j)
            
    return distanceE, saveIndex


def fittingBest(x,y):
    #variables 
    distanceE = []
    saveIndex = []
    
    # determine best fit line
    #TODO: check issuhere
    par = np.polyfit(x,y, 1, full=True)
    
    #use np.polyfit
    #https://stackoverflow.com/questions/18767523/fitting-data-with-numpy/18767992#18767992
    
    
def reOrderPoints(x,y):
    k = list(zip(y,x)) #backwards order since this is copied from roughness script
    
    #find starting point of contour
    minIndices = np.where(x == k[:,1].min())[0]
    minPoints = k[minIndices]
    minIndx = np.where(minPoints[:,0] == minPoints[:,0].min())[0][0]
    startingCord = k[minIndices[minIndx]]
    
    #array to store ordered points
    newOrder = [startingCord]
    
    #delete starting point from contour array (only pairs values in k)
    k = np.delete(k, minIndx, axis=0)
    
    
    #Find nearest neighbour, stop when next vertex is dist > 4 away
    while(len(k) > 1):
        nbrs = NearestNeighbors(n_neighbors=2, algorithm='ball_tree').fit(k)
        distance, indices = nbrs.kneighbors([newOrder[-1]])
        
        if(distance[0][0] > 4):
            break
        else:
            indices = indices[:,0]
            newOrder.append(k[indices[0]])
            k = np.delete(k, indices[0], axis=0)


    #get unqiue points, maintain order
    _, idx = np.unique(newOrder, axis=0,  return_index=True)
    newOrderIndx = np.sort(idx)
    
    finalOrder = []
    
    for p in newOrderIndx:
        finalOrder.append(newOrder[p])
        
    
    finalOrder = np.array(finalOrder)
    
    return  np.array(finalOrder[:,0]), np.array(finalOrder[:,1])

#===============GLOBAL VARIABLES=================================
inputDir = "C:/Users/v.jayaweera/Documents/Anne/Line-Fitting/LineProfile-Output/Probe3"

Ra = []


for file in os.listdir(inputDir):
    degree = file.split(".")[0].split("_")[2]
    
    distanceE =[]
    saveIndex = []
    
    data = np.loadtxt(inputDir + "/"+ file)
    r, z = data[:,0], data[:,1]
    
    
    #create kernel
    sig = 350
    size = 319
    kernel = fb.gauss1D(size, sig)   

    #get baseline - GAUSS
    rscipy = signal.convolve(r, kernel, mode='valid')
    zscipy = signal.convolve(z, kernel, mode='valid')
    
    dx = np.diff(rscipy)
    dy = np.diff(zscipy)
    
    
    #Plotting all data
    plt.plot(rscipy,zscipy, "r.-")
    plt.plot(r,z, "b.-")
    plt.show()
  
    finalOrder = list(zip(r,z))
    polyGon = shapely.geometry.LineString(finalOrder)
     
    for j in range(1,len(dx)):
        xs, ys = fb.createNormalLine(rscipy[j], zscipy[j], dx[j], dy[j])
       
        
        stack = np.stack((xs,ys), axis=-1)
        line = shapely.geometry.LineString(stack)
        
        #TODO remove this from main CODE
        if(polyGon.intersects(line)):
            #intersection geometry
            interPoints = polyGon.intersection(line)
            
            #intersection point
            mx, my = fb.proccessIntersectionPoint(interPoints, rscipy[j], zscipy[j])
            
            euD = fb.euclidDist(rscipy[j], zscipy[j], mx, my)
            distanceE.append(euD)
            saveIndex.append(j)
    
    Ra.append([distanceE[-1], int(degree)])

    # # calculations 
    # mean = np.average(distanceE)
    
   
Ra = np.array(Ra)

#save to file 
df = pd.DataFrame(data=Ra, columns=['Ra', 'Degree' ])
df.to_csv("Probe3_Gauss_Ra_Data.csv", index=False)
    