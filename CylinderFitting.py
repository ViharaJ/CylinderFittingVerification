"""
POINT CLOUD name syntax: Side_ProbeX_*.txt where X is an int

Fits cylinder to point cloud
Writes r,phi files to folder name LineProfile-Output 
    - this folder must be created in the directory of this script 
    
    
Credit to answer: https://stackoverflow.com/questions/43784618/fit-a-cylinder-to-scattered-3d-xyz-point-data
"""
import numpy as np
import os
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
import csv
import pandas as pd


#===========================FUNCTIONS=================================
def rotX(data, angle):
    '''
    angle: in radians
    '''
    rotXmatrix = np.array([[1,0,0],
                           [0, np.cos(angle), -np.sin(angle)],
                           [0, np.sin(angle), np.cos(angle)]])
    

    
    data_rot = np.dot(rotXmatrix, np.transpose(data))
    
    return np.transpose(data_rot)


def rotY(data, angle):
    '''
    angle: in radians
    '''
    rotYmatrix = np.array([[np.cos(angle), 0, np.sin(angle)],
                           [0, 1, 0],
                           [-np.sin(angle), 0, np.cos(angle)]])
    
    data_rot = np.dot(rotYmatrix, np.transpose(data))
    
    return np.transpose(data_rot)


def cart2pol(x, y):
    '''
    returns r and phi, angle is in degrees [-180,180]
    '''
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)* 180 / np.pi
    return(rho, phi)


def translate(data, x, y):
    '''
    Translate data by x and y to (0,0)
    '''
    data[:,0] = data[:,0] - x
    data[:,1] = data[:,1] - y
    
    return data


def writeToTxt(name, a, b):
    d  = np.column_stack([a, b])
    np.savetxt('LineProfile-Output/' + name + '.txt', d, fmt=['%f','%f'])
    
    
def cylinderFitting(xyz,p,th):
    """
    This is a fitting for a vertical cylinder fitting
    Reference:
    http://www.int-arch-photogramm-remote-sens-spatial-inf-sci.net/XXXIX-B5/169/2012/isprsarchives-XXXIX-B5-169-2012.pdf

    xyz is a matrix contain at least 5 rows, and each row stores x y z of a cylindrical surface
    p is initial values of the parameter;
    p[0] = Xc, x coordinate of the cylinder centre
    P[1] = Yc, y coordinate of the cylinder centre
    P[2] = alpha, rotation angle (radian) about the x-axis
    P[3] = beta, rotation angle (radian) about the y-axis
    P[4] = r, radius of the cylinder

    th, threshold for the convergence of the least squares

    """   
    x = xyz[:,0]
    y = xyz[:,1]
    z = xyz[:,2]

    fitfunc = lambda p, x, y, z: (- np.cos(p[3])*(p[0] - x) - z*np.cos(p[2])*np.sin(p[3]) - np.sin(p[2])*np.sin(p[3])*(p[1] - y))**2 + (z*np.sin(p[2]) - np.cos(p[2])*(p[1] - y))**2 #fit function
    errfunc = lambda p, x, y, z: fitfunc(p, x, y, z) - p[4]**2 #error function 

    est_p , success = leastsq(errfunc, p, args=(x, y, z), maxfev=1000)

    return est_p


def data_for_cylinder_along_z(center_x,center_y,radius,height_z):
    z = np.linspace(0, height_z, 50)
    theta = np.linspace(0, 2*np.pi, 50)
    theta_grid, z_grid=np.meshgrid(theta, z)
    x_grid = radius*np.cos(theta_grid) + center_x
    y_grid = radius*np.sin(theta_grid) + center_y
    return x_grid,y_grid,z_grid


def deg2rad(angle):
    return angle *( np.pi / 180)

# PLOTTING STUFF

# CT data and fitted cylinder
# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.scatter(relv_x, relv_y, relv_z, marker=".", color="green")
# ax.plot_surface(Xc, Yc, Zc, alpha=0.8) #plot fitted cylinder
# ax.set_title('3D line plot geeks for geeks')
# plt.show()

# R vs Z
# plt.plot(relv_z, relv_radii, '.')
# plt.xlabel("z")
# plt.ylim([0,1])
# plt.ylabel("r")
# plt.show()

# adjusting axes
# ax.axes.set_xlim3d (left=-3, right=2) 
# ax.axes.set_ylim3d (bottom=-3, top=2) 
# ax.axes.set_zlim3d (bottom=0, top=3) 

# PLOT FITTED CYLINDER, center = 0,0
# X0.c,Yc,Zc = data_for_cylinder_along_z(0,0,est_p[4], 5)
# ax.plot_surface(Xc, Yc, Zc, alpha=8)
#===============GLOBAL VARIABLES=================================
inputDir = "C:/Users/v.jayaweera/Documents/Anne/Line-Fitting/Point-Cloud-Files"
stepWidth = 1

#=======================MAIN===================================
for file in os.listdir(inputDir):
    probe = file.split("_")[1]
    #include if statement to filter inocrrect file types here 
    data = np.loadtxt(inputDir + "/"+ file, skiprows=1)
    
    #shift to z = 0
    data[:,2] = data[:,2] - ((np.max(data[:, 2]) + np.min(data[:,2]))/2)
    
    # fit cylinder
    p = np.array([0,0,0,0,0.8]) # initial fit parametrs
    est_p =  cylinderFitting(data,p,0.000001)
    print ("Fitting Done!\n")
    print ("Estimated Parameters for ", file, ": ")
    print (est_p)
    
    #Translate CT data to x,y,z = 0
    data = translate(data, est_p[0], est_p[1])
    
    print("Angles ", est_p[2]*( 180 / np.pi), est_p[3]* (180 / np.pi))
   
    # #ROTATION
    # data = rotX(data, -(est_p[2]))
    # data = rotY(data, deg2rad(90))
    
    x = data[:,0][::100]
    y = data[:,1][::100]
    z = data[:,2][::100]
    
    # #convert to cylindrical coordinates
    # (r, theta) = cart2pol(data[:,0], data[:,1])
    
    # theta = theta + 180 #change range to [0, 360]
    
    # # PLOTTING
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    
    
    # # PLOT ORIGINAL
    ax.scatter(x, y, z, marker=".", color="green")    
    ax.set_title(file)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    
    
    #Plot a fitted cylinder 
    Xc,Yc,Zc = data_for_cylinder_along_z(0,0,est_p[4], 5)
    fit_data = np.transpose(np.array([Xc.flatten(), Yc.flatten(), Zc.flatten()]))                           
  
    
    fit_data[:,2] = fit_data[:,2] - ((np.max(fit_data[:, 2]) + np.min(fit_data[:,2]))/2)
    fit_data = rotX(fit_data, est_p[2])
    fit_data = rotX(fit_data, est_p[3])

    ax.scatter(fit_data[:,0],fit_data[:,1],fit_data[:,2], marker=".", color="red") 
    plt.show()
        
        
    
   




#====================TESTING=-========================

# mat1 = np.array([[1,2,3], [4,5,6], [6,7,8], [6,7,8]])
 
mat2 = [[1,2,5], [2,2,6], [3,2,7], [4,2,8], [5,2,3]]

# print(rotX(mat1, 0.52359878))

# print(rotX(mat1, 0.95993109))