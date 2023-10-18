"""
Fits cylinder to point cloud
Credit to answer: https://stackoverflow.com/questions/43784618/fit-a-cylinder-to-scattered-3d-xyz-point-datahttps://stackoverflow.com/questions/43784618/fit-a-cylinder-to-scattered-3d-xyz-point-data
"""
import numpy as np
import os
from scipy.optimize import leastsq
import matplotlib.pyplot as plt


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

def cart2pol(x, y):
    '''
    returns r and phi, angle is in degrees [-180,180]
    '''
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)* 180 / np.pi
    return(rho, phi)


#===============GLOBAL VARIABLES=================================
inputDir = "C:/Users/v.jayaweera/Documents/Anne/Line-Fitting/Point-Cloud-Files"
stepWidth = 1


#=======================MAIN===================================
for file in os.listdir(inputDir):
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
    x = data[:,0][::150]
    y = data[:,1][::150]
    z = data[:,2][::150]
    
    #convert to cylindrical coordinates
    (r, theta) = cart2pol(data[:,0], data[:,1])
    
    # PLOTTING
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    
    # PLOT ORIGINAL
    ax.scatter(x, y, z, marker=".", color="green")
    
    # PLOT FITTED CYLINDER
    Xc,Yc,Zc = data_for_cylinder_along_z(est_p[0], est_p[1],est_p[4], 5)
    ax.plot_surface(Xc, Yc, Zc, alpha=0.8)
    
    ax.set_title('3D line plot geeks for geeks')
    plt.show()
    
   
    
    for deg in range(-180,180, 1):
        indices = np.where(np.logical_and(theta > deg - stepWidth/2, theta <= deg + stepWidth/2))

        relv_x = data[indices, 0]
        relv_y = data[indices, 1]
        relv_z = data[indices, 2]
        
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.scatter(relv_x, relv_y, relv_z, marker=".", color="green")
        ax.set_title('3D line plot geeks for geeks')
        plt.show()