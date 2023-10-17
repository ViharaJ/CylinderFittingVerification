"""
Created on Thu Oct 12 15:13:35 2023

Get points of line

Create baseline using 2 methods 
    1. Linear regression 
    2. Gaussian kernel 



"""
import os 
import cv2 
import numpy as np
import Module.Functions as fb
import shapely
from scipy import signal

inputDir = ""

acceptedFileType = None


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
        
    

for file in os.listdir(inputDir):
    #include if statement to filter inocrrect file types here 
    
    #TODO: Extract from point cloud
    x,y = []
    
    
    #create kernel
    sig = 350
    size = 319
    kernel = fb.gauss1D(size, sig)   

    #get baseline - GAUSS
    xscipy = signal.convolve(x, kernel, mode='valid')
    yscipy = signal.convolve(y, kernel, mode='valid')
    
    dx = np.diff(xscipy)
    dy = np.diff(yscipy)

    #best fit 
    # determine best fit line
    par = np.polyfit(x,y, 1, full=True)
  
    
    
    # #TODO: place best fit line points
    # bestFitLine = shapely.geometry.LineString(list(zip(x,y)))
    
    
    # for i in range(len(x)):
    #     xs, ys = fb.createNormalLine(x,y, dx[i], dy[i])
    #     stack = np.stack((xs,ys), axis=-1)
    #     line = shapely.geometry.LineString(stack)
        
    #     #TODO remove this from main CODE
    #     if(bestFitLine.intersects(line)):
    #         #intersection geometry
    #         interPoints = bestFitLine.intersection(line)
            
    #         #intersection point
    #         mx, my = fb.proccessIntersectionPoint(interPoints, xscipy[j], yscipy[j])
    