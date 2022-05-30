import matplotlib.pyplot as plt
import numpy as np
import random
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay, ConvexHull
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import math

def deltapeak(x,sigma):
    delta = np.exp(-(x)**2/(2*sigma))/np.sqrt(2*np.pi*sigma)
    return delta


def brdf_lambertian(rhod,theta_i,theta_r,phi_i,phi_r):
    f = np.zeros(2*res*res)
    for i in range(2*res*res):
            f[i] = rhod/np.pi
    return f



def brdf_specular3(rhod,theta_i,theta_r,phi_i,phi_r,R_i,sigma):
    z = rhod*R_i*(deltapeak(theta_r-theta_i,sigma)*deltapeak(phi_r-phi_i-np.pi,sigma)/(np.cos(theta_r)*np.sin(theta_r)+epsilon))
    return z

def brdf_specular4(rhod,theta_i,theta_r,phi_i,phi_r,R_i,sigma):
    z = rhod*R_i*(deltapeak(np.cos(theta_r)-np.cos(theta_i),sigma)*deltapeak(phi_r-phi_i-np.pi,sigma)/(np.cos(theta_r)+epsilon))
    return z


def brdf_lambertian2(rhod,theta_i,theta_r,phi_i,phi_r):
    z = rhod*(np.cos(theta_i)*np.cos(theta_r))**0/np.pi
    return z
def brdf_opposition(rhod,theta_i,theta_r,k):
    z =  rhod*(np.cos(theta_i)*np.cos(theta_r))**k/np.pi
    return z

def brdf_opposition2(rhod,theta_i,theta_r,k):
    z =  rhod*(np.cos(theta_i)*np.cos(theta_r)+epsilon)**(-k)/np.pi
    return z



epsilon = 0.001
sigma = 0.01
R_i = 0.8
rhod = 0.8
res = 50
k = 0.5

theta_i = np.pi/3
phi_i = np.pi/4

phi_r = np.linspace(0,2*np.pi,2*res)
theta_r = np.linspace(0,np.pi/2,res)


theta_i_opp = np.linspace(0,np.pi/2,res)

phi_mesh, theta_mesh = np.meshgrid(phi_r,theta_r)

theta_imeshopp, theta_rmeshopp = np.meshgrid(theta_i_opp, theta_r)



brdfLamb = brdf_lambertian(rhod,theta_i,theta_r,phi_i,phi_r)
brdfLamb_grid = brdfLamb.reshape((2*res,res),order = 'F')

brdfLamb2 = brdf_lambertian2(rhod,theta_i,theta_mesh,phi_i,phi_mesh)


brdfSpec3 = brdf_specular3(rhod,theta_i,theta_mesh,phi_i,phi_mesh,R_i,sigma)
brdfSpec4 = brdf_specular4(rhod,theta_i,theta_mesh,phi_i,phi_mesh,R_i,sigma)


brdfopp1 = brdf_opposition(rhod,theta_imeshopp,theta_rmeshopp,k)
brdfopp2 = brdf_opposition2(rhod,theta_imeshopp,theta_rmeshopp,k)



## plotting the brdf's

##fig = plt.figure(figsize=(13, 7))
##ax = plt.axes(projection='3d')
##surf = ax.plot_surface(phi_mesh,theta_mesh, brdfLamb2, rstride=1, cstride=1, cmap='gray', edgecolor='none')
##ax.set_xlabel('phi_r')
##ax.set_ylabel('theta_r')
##ax.set_zlabel('brdf')
##ax.set_title('Lambertian, phi_i = pi/4, theta_i = pi/3')
##fig.colorbar(surf, shrink=0.7, aspect=10) 
##ax.view_init(90, 0)
##             



##fig = plt.figure(figsize=(13, 7))
##ax = plt.axes(projection='3d')
##surf = ax.plot_surface(phi_mesh,theta_mesh, brdfSpec3, rstride=1, cstride=1, cmap='gray', edgecolor='none')
##ax.set_xlabel('phi_r')
##ax.set_ylabel('theta_r')
##ax.set_zlabel('brdf')
##ax.set_title('Specular2, theta_i = pi/2, phi_i = pi/2')
##fig.colorbar(surf, shrink=0.7, aspect=10) 
##ax.view_init(90, 0)

fig = plt.figure(figsize=(13, 7))
ax = plt.axes(projection='3d')
surf = ax.plot_surface(phi_mesh,theta_mesh, brdfSpec4, rstride=1, cstride=1, cmap='gray', edgecolor='none')
ax.set_xlabel('phi_r')
ax.set_ylabel('theta_r')
ax.set_zlabel('brdf')
ax.set_title('Specular2, phi_i = pi/4, theta_i = pi/3, sigma='+str(sigma))
fig.colorbar(surf, shrink=0.7, aspect=10) 
ax.view_init(90, 0)

##fig = plt.figure(figsize=(13, 7))
##ax = plt.axes(projection='3d')
##surf = ax.plot_surface(theta_imeshopp,theta_rmeshopp, brdfopp1, rstride=1, cstride=1, cmap='gray', edgecolor='none')
##ax.set_xlabel('theta_i')
##ax.set_ylabel('theta_r')
##ax.set_zlabel('brdf')
##ax.set_title('Opposition surge, k = '+str(k))
##fig.colorbar(surf, shrink=0.7, aspect=10) 
##ax.view_init(90, 0)

##fig = plt.figure(figsize=(13, 7))
##ax = plt.axes(projection='3d')
##surf = ax.plot_surface(theta_imeshopp,theta_rmeshopp, brdfopp2, rstride=1, cstride=1, cmap='gray', edgecolor='none')
##ax.set_xlabel('theta_i')
##ax.set_ylabel('theta_r')
##ax.set_zlabel('brdf')
##ax.set_title('New Opposition surge, k = -'+str(k))
##fig.colorbar(surf, shrink=0.7, aspect=10) 
##ax.view_init(90, 0)


plt.show()
