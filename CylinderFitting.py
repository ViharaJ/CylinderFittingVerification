"""
Fits cylinder to point cloud
Writes r,phi files to folder name LineProfile-Output 
    - this folder must be created in the directory of this script 
Credit to answer: https://stackoverflow.com/questions/43784618/fit-a-cylinder-to-scattered-3d-xyz-point-datahttps://stackoverflow.com/questions/43784618/fit-a-cylinder-to-scattered-3d-xyz-point-data
"""
import numpy as np
import os
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
import csv
import pandas as pd


#===========================FUNCTIONS=================================
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


def translate(data, x, y):
    '''
    Translate data by x and y
    Also shift to z = 0
    '''
    data[:,0] = data[:,0] - x
    data[:,1] = data[:,1] - y
    data[:,2] = data[:,2] - np.min(data[:,2])
    
    return data

def rotX(data, angle):
    rotXmatrix = np.array([[1,0,0],
                           [0, np.cos(angle), -np.sin(angle)],
                           [0, np.sin(angle), np.cos(angle)]])
    
    return np.dot(data, rotXmatrix)


def rotY(data, angle):
    rotYmatrix = np.array([[np.cos(angle), 0, np.sin(angle)],
                           [0, 1, 0],
                           [-np.sin(angle), 0, np.cos(angle)]])
    
    return np.dot(data, rotYmatrix)

def cart2pol(x, y):
    '''
    returns r and phi, angle is in degrees [-180,180]
    '''
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)* 180 / np.pi
    return(rho, phi)



def writeToTxt(name, a, b):
    d  = np.column_stack([a, b])
    np.savetxt('LineProfile-Output/' + name + '.txt', d, fmt=['%f','%f'])
    
        

#===============GLOBAL VARIABLES=================================
inputDir = "C:/Users/v.jayaweera/Documents/Anne/Line-Fitting/Point-Cloud-Files"
stepWidth = 1


#=======================MAIN===================================
for file in os.listdir(inputDir):
    probe = file.split("_")[1]
    #include if statement to filter inocrrect file types here 
    data = np.loadtxt(inputDir + "/"+ file, skiprows=1)
    
    # fit cylinder
    p = np.array([-13.79,-8.45,0,0,0.3])
    est_p =  cylinderFitting(data,p,0.00001)
    print ("Fitting Done!\n")
    print ("Estimated Parameters: ")
    print (est_p)
    
    #shift CT data to x,y,z = 0
    data = translate(data, est_p[0], est_p[1])
    #ROTATION
    data = rotX(data, est_p[2])
    data = rotY(data, est_p[3])
    
    x = data[:,0][::150]
    y = data[:,1][::150]
    z = data[:,2][::150]
    
    #convert to cylindrical coordinates
    (r, theta) = cart2pol(data[:,0], data[:,1])
    
    theta = theta + 180 #change range to [0, 360]
    # PLOTTING
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    # ax.axes.set_xlim3d (left=-3, right=2) 
    # ax.axes.set_ylim3d (bottom=-3, top=2) 
    # ax.axes.set_zlim3d (bottom=0, top=3) 
    # PLOT ORIGINAL
    ax.scatter(x, y, z, marker=".", color="green")
    
    # PLOT FITTED CYLINDER, center = 0,0
    Xc,Yc,Zc = data_for_cylinder_along_z(0,0,est_p[4], 5)
    ax.plot_surface(Xc, Yc, Zc, alpha=0.8)
    
    ax.set_title('3D Plot')
    plt.show()
    
   
    Ra = []
    for deg in range(0,360, 1):
        indices = np.where(np.logical_and(theta > deg - stepWidth/2, theta <= deg + stepWidth/2))

        relv_x = data[indices, 0]
        relv_y = data[indices, 1]
        relv_z = data[indices, 2][0]
        relv_radii = r[indices]
        
        # plt.plot(relv_z, relv_radii, '.')
        # plt.xlabel("z")
        # plt.ylim([0,1])
        # plt.ylabel("r")
        # plt.show()
        
        writeToTxt(probe + "/Profile_Deg_" + str(deg), relv_z, relv_radii)
        
        # plot splice
        # fig = plt.figure()
        # ax = plt.axes(projection='3d')
        # ax.scatter(relv_x, relv_y, relv_z, marker=".", color="green")
        # ax.plot_surface(Xc, Yc, Zc, alpha=0.8) #plot fitted cylinder
        # ax.set_title('3D line plot geeks for geeks')
        # plt.show()
        
        Ra.append([np.mean(np.abs(relv_radii - np.mean(relv_radii))), deg])
        
    Ra = np.array(Ra)
    
    #save to file 
    df = pd.DataFrame(data=Ra, columns=['Ra', 'Degree' ])
    df.to_csv(probe + "_Ra_Data.csv", index=False)

    
    plt.polar(Ra[:,1] / (180 / np.pi), Ra[:,0],'.')
    plt.show()
    
    plt.plot(Ra[:,1], Ra[:,0], 'b.-')
    plt.xlabel("Angel")
    plt.ylabel("rho")
    plt.show()