# -*- coding: utf-8 -*-
"""
Function reorders points
"""
import numpy as np
from sklearn.neighbors import NearestNeighbors

def reOrderPoints(x,y, rem_duplicates = False):
    '''
    EXPECTED INPUT: profile is orientated horizontally
    returns: ordered points if sucessful, None if failed
    '''
    #turn to array of shape (n,2)
    k = np.array(list(zip(x,y)))
    
    #length for comparison of reordering
    len_comparison = (len(k)/2)*0.95 if rem_duplicates else len(k)*0.95
    

    #find starting point of contour
    minIndices = np.where(k[:,1] == k[:,1].min())[0]
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
    if rem_duplicates:
        _, idx = np.unique(newOrder, axis=0,  return_index=True)
        newOrderIndx = np.sort(idx)
    
        finalOrder = []
        
        for p in newOrderIndx:
            finalOrder.append(newOrder[p])
    else: 
        finalOrder = newOrder
        
    
    finalOrder = np.array(finalOrder)
    
    if(len(finalOrder) >= len_comparison):                    
        x = np.array(finalOrder[:,0])
        y = np.array(finalOrder[:,1])
        
        return x,y
    
    return None
    
