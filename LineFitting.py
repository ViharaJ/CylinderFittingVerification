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

#===============GLOBAL VARIABLES=================================
inputDir = "C:/Users/v.jayaweera/Documents/Anne/Line-Fitting/Point-Cloud-Files"

acceptedFileType = None


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
    



for file in os.listdir(inputDir):
    