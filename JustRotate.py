import os
import numpy as np
import matplotlib.pyplot as plt

inputDir = "C:/Users/v.jayaweera/Documents/Anne/Line-Fitting/Point-Cloud-Files/Side_Probe3_Cylinder1_ObjCoord.txt"
stepWidth = 1


def rotX(data, angle):
    '''
    angle: in radians
    '''
    rotXmatrix = np.array([[1,0,0],
                           [0, np.cos(angle), -np.sin(angle)],
                           [0, np.sin(angle), np.cos(angle)]])
    

    
    data_rot = np.dot(rotXmatrix, np.transpose(data))
    
    return np.transpose(data_rot)


def deg2rad(angle):
    return angle *( np.pi / 180)

#=======================MAIN===================================

#include if statement to filter inocrrect file types here 
data = np.loadtxt(inputDir, skiprows=1)

#shift to z = 0
data[:,2] = data[:,2] - ((np.max(data[:, 2]) + np.min(data[:,2]))/2)
data[:,1] = data[:,1] - ((np.max(data[:, 1]) + np.min(data[:,1]))/2)
data[:,0] = data[:,0] - ((np.max(data[:, 0]) + np.min(data[:,0]))/2)

data = rotX(data, deg2rad(-90))

x = data[:,0][::50]
y = data[:,1][::50]
z = data[:,2][::50]


# # PLOTTING
fig = plt.figure()
ax = plt.axes(projection='3d')

ax.set_aspect('equal')
# # PLOT ORIGINAL
ax.scatter(x, y, z, marker=".", color="green")    

ax.set_xlabel("x")
ax.set_ylabel("y")

plt.show()