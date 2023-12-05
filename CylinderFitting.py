"""
Reads point cloud file and finds the roughness along vertical profiles. 
The script creates an excel file per file, indicating the roughness at each angle

2 Directories are created inside the directory containing the point cloud files.
    1. LineProfile-Output contains the excel documents
    2. Figures contains the plots


How to use:
    1. Place script in it's own directory
    2. Change inputDir to location of Point Cloud files 
    

How it works: 
    1. The point cloud file is loaded 
    2. An ideal cyclinder is fitted to it
        Credit to answer: 
            https://stackoverflow.com/questions/43784618/fit-a-cylinder-to-scattered-3d-xyz-point-data
    3. Using the parameters of the ideal cyclinder, reorient the cyclinder 
        to be upright and with no rotations 
    4. Iterate over the angles in the range [0,360]
    5. At each angle, extract a vertical profile. 
    6. Do a linear fit the data 
    7. Compute roughness between linear fitted line and surface
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
    converts cartesian coordinates into polar coordinates
    Abstand vom Pol - r / rho (Radius, Radialkoordinate)
    Winkel -  phi (Winkelkoordinate)
  
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

        
def cylinderFitting(xyz,p,th):
    """
    This is a fitting for a vertical cylinder fitting
    Reference:
    http://www.int-arch-photogramm-remote-sens-spatial-inf-sci.net/XXXIX-B5/169/2012/isprsarchives-XXXIX-B5-169-2012.pdf
    https://stackoverflow.com/questions/43784618/fit-a-cylinder-to-scattered-3d-xyz-point-data/44164662#44164662

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

    fitfunc = lambda p, x, y, z: (- np.cos(p[3])*(p[0] - x) - z*np.cos(p[2])*np.sin(p[3]) - 
                                  np.sin(p[2])*np.sin(p[3])*(p[1] - y))**2 + (z*np.sin(p[2]) - np.cos(p[2])*(p[1] - y))**2 #fit function
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


def createDir(root, folderName): 
    '''
    creates new folder if it doesn't exist
    returns: new folder's path
    '''
    newPath = os.path.join(root, folderName)
    if not os.path.exists(newPath):
        os.makedirs(newPath)
        
    return newPath


#=============================GLOBAL VARIABLES=================================
inputDir = "Z:/Projekte/42029-FOR5250/Vihara/Line-Fitting/Point-Cloud-Files" #CHANGE TO YOUR DIRECTORY
stepWidth = 1


# create directories to save files
savePath = createDir(inputDir, "LineProfile-Output")
figDir = createDir(inputDir, "Figures")


#===================================MAIN===================================
for file in os.listdir(inputDir):
    #include if statement to filter incorrect file types here 
    if file.split(".")[-1] == "txt":
        data = np.loadtxt(inputDir + "/"+ file, skiprows=1)
        
        Ra = []
        angle = []
        
        #shift center of point cloud to x, y, z = 0
        data[:,0] = data[:,0] - ((np.max(data[:, 0]) + np.min(data[:,0]))/2)
        data[:,1] = data[:,1] - ((np.max(data[:, 1]) + np.min(data[:,1]))/2)
        data[:,2] = data[:,2] - ((np.max(data[:, 2]) + np.min(data[:,2]))/2)
        
        # fit cylinder    
        p = np.array([0,0,0,0,0.8]) # initial fit parametrs
        '''
        p[0] = Xc, x coordinate of the cylinder centre
        P[1] = Yc, y coordinate of the cylinder centre
        P[2] = alpha, rotation angle (radian) about the x-axis
        P[3] = beta, rotation angle (radian) about the y-axis
        P[4] = r, radius of the cylinder - initial value based on design 
        '''
        est_p =  cylinderFitting(data,p,0.00001)
        print ("Fitting Done!\n")
        print ("Estimated Parameters for ", file, ": ")
        print (est_p)
        
        #Translate CT data to x,y = 0
        data = translate(data, est_p[0], est_p[1])
        
        # in radians and converting to degrees 
        print("Angles ", est_p[2]*( 180 / np.pi), est_p[3]* (180 / np.pi))
       
        # Rotate to be upright, 
        # TODO: check if suitable
        data = rotX(data, -(est_p[2]))
        data = rotY(data, -(est_p[3]))
        
        x = data[:,0]
        y = data[:,1]
        z = data[:,2]
        
        #convert to cylindrical coordinates
        # Abstand vom Pol - r (Radius), theta - Winkelkoordinate
        (r, theta) = cart2pol(data[:,0], data[:,1])
        
        
        # Create a fitted cylinder > used for the 3D plot not for surface calculation
        Xc,Yc,Zc = data_for_cylinder_along_z(0,0,est_p[4], 2)
        fit_data = np.transpose(np.array([Xc.flatten(), Yc.flatten(), Zc.flatten()]))  
    
        
        # iterate over the angles (different name than theta, because phi is the loop variable and would otherwise overwrite theta)
        for phi in range(-180, 180, stepWidth): #-180 to 180 because cart2pol defines it in this range
        
        
            #get indices where theta is phi
            indices = np.where(np.logical_and(theta >= (phi - stepWidth/2), 
                                                        theta < (phi+stepWidth / 2)))
            interval_r = r[indices]
            relv_x = x[indices]  
            relv_y = y[indices] 
            relv_z = z[indices] # used for printing, not calculation
            
            # find vals for line of best fit
            a, b = np.polyfit(relv_z, interval_r, 1)
            
            # line of best fit 
            expected_rho = a*relv_z+b
            
            # roughness 
            Ra.append(np.mean(np.abs(interval_r - expected_rho)) * 1)
            angle.append(phi)
        
        
         
        # plot Ra vs angle data 
        plt.clf()
        plt.plot(angle, Ra, "g.")
        plt.show()
        
        # Plot 3D figures
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        
        # Plot origina point cloud 
        ax.scatter(x, y, z, marker=".", color="green")    
        ax.set_title(file)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        
    
        # plot ideal cyclinder 
        # TODO: DEBUG WHY IT LOOKS WEIRD
        fit_data[:,2] = fit_data[:,2] - ((np.max(fit_data[:, 2]) + np.min(fit_data[:,2]))/2)
        ax.scatter(fit_data[:,0],fit_data[:,1],fit_data[:,2], marker=".", color="red") 
        plt.show()
        
        
        # save data
        d = {'Angle': angle, 'Ra(microns)': Ra}
        df = pd.DataFrame(data=d)
        df.to_excel(os.path.join(savePath, file.split('.')[0] +  '.xlsx'), index=False)
        
        # plot polar plot and save
        plt.axes(projection = 'polar') 
        plt.polar(deg2rad(np.array(angle)), Ra) 
        plt.savefig(os.path.join(figDir, file.split(".")[0]), dpi=1000)
        plt.show()
            
        #THE PORTION BELOW CAN BE COMMENTED OUT
        # Plot levels, cross sections
        vals = np.linspace(-0.7, 0.5, 100)
        
        for i in range(len(vals)):
            level = vals[i]
            indices = np.where(np.logical_and(z > level-0.01, z < level+0.01))
            
            plt.title(file)
            plt.plot(x[indices], y[indices], "g.")
            plt.show()
        
   
    
    

