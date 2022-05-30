# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import random
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay, ConvexHull
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap



def create_sphere1(r, res):
    '''
    creates points on a sphere in the form:
        [[ x1.  y1.  z1.]
         ....
         [ xn.  yn.  zn.]]
    '''
    phi = np.linspace(0, 2*np.pi, 2*res)
    theta = np.linspace(0, np.pi, res)

    theta, phi = np.meshgrid(theta, phi)  
    
    theta_flat = np.ndarray.flatten(theta, order = 'F')
    phi_flat = np.ndarray.flatten(phi, order = 'F')
    
    r_pre = r*np.sin(theta_flat)
    x = np.cos(phi_flat)*r_pre
    y = np.sin(phi_flat)*r_pre
    z = r*np.cos(theta_flat)

    ''' 
    flatten arrays, delete top- end bottom points of circle and append both points once
    '''
#    x = np.delete(x, [0, res-1], 1)
#    x = np.ndarray.flatten(x, order ='F') 
#    x = np.append(x, [0,0])
#    y = np.delete(y, [0, res-1], 1)
#    y = np.ndarray.flatten(y, order = 'F') 
#    y= np.append(y, [0,0])
#    z = np.delete(z, [0, res-1], 1)
#    z = np.ndarray.flatten(z, order = 'F') 
#    z= np.append(z, [1,-1])
    

    coordinates = np.vstack((x,y,z)).transpose()

    return coordinates, theta, phi





def pnt_in_cvex_hull(hull, pnt):
    '''
    Checks if `pnt` is inside the convex hull.
    '''
    new_hull = ConvexHull(np.concatenate((hull.points, [pnt])))
    if np.array_equal(new_hull.vertices, hull.vertices): 
        return True
    else:
        return False




def height(p, vertices, seeds, altitudes):
    max_dist = 0
    '''
    finding shortest edge
    '''
    for i in range(3):
        for j in range(i+1 , 4):
            dist = np.linalg.norm(vertices[i] - vertices[j])
            if dist > max_dist:
                max_dist = dist
                furthest = np.array([i, j])   
   
    '''
    reordering everything
    '''          
    vertices[[0, furthest[0]]] = vertices[[furthest[0], 0]]
    vertices[[1, furthest[1]]] = vertices[[furthest[1], 1]]
    
    seeds[[0,furthest[0]]] = seeds[[furthest[0], 0]]
    seeds[[1,furthest[1]]] = seeds[[furthest[1], 1]]
    
    altitudes[[0, furthest[0]]] = altitudes[[furthest[0], 0]]
    altitudes[[1, furthest[1]]] = altitudes[[furthest[1], 1]]
    
    '''
    creating new edge
    new altitude = average altitude + random(-0,05; 0,05)*sqrt(distance)
    '''
    v_new = (vertices[0]+vertices[1])/2
    s_new = (seeds[0]+seeds[1])/2
    random.seed(s_new)
    a_new = (altitudes[0]+ altitudes[1])/2 + 0.01*(random.random()- 0.5)*max_dist**1.5

    
    '''
    finding in which tetrahedon our point p is
    '''    
    tetra = np.copy(vertices)
    tetra[1] = v_new
    hull = ConvexHull(tetra)
    
    if pnt_in_cvex_hull(hull, p):
        vertices[1] = v_new
        seeds[1] = s_new
        altitudes[1] = a_new
    else:
        vertices[0] = v_new
        seeds[0] = s_new
        altitudes[0] = a_new 
    
    '''
    stop if resolution is great enough
    '''
    if max_dist < 0.001:
#        print(np.sum(altitudes)/4)
        return np.sum(altitudes)/4
    else:
#        print(vertices)
#        print(seeds)
#        print(altitudes)
        
        height(p, vertices, seeds, altitudes)
    return np.sum(altitudes)/4




'''
initialize sphere and final altitudes vector
'''
radius = 1
res = 10
sphere, theta, phi = create_sphere1(radius,res)
alts = np.zeros(len(sphere[:, 0]))
albedo_map = np.zeros(len(sphere[:, 0]))
'''
initialize vertices, seeds and altitudes(sea level)
'''


'''
calculate altitude for every point on sphere
'''
for i in range(len(sphere[:, 0])):
    vertices_ = np.array([[-2  ,  2.1,  2.2], [ 2  , -2.1,  2.2], [ 2  ,  2.1, -2.2],[-2.3,  -2.4,  -2.5]])
    seeds_ , altitudes_ = np.array([30,42,1,500]), np.array([0., 0. , 0., 0.])
    point = sphere[i, :]
    alt = height(point, vertices_, seeds_, altitudes_)
    alts[i] = alt
#    print(i, ': ', alt)

print('busy')
'''
create corresponding meshgrid
'''
alt_grid = alts.reshape((2*res, res), order = 'F')

#
fig = plt.figure(figsize=(13, 7))
ax = plt.axes(projection='3d')
surf = ax.plot_surface(theta,phi, alt_grid, rstride=1, cstride=1, cmap='gray', edgecolor='none')
ax.set_xlabel('theta')
ax.set_ylabel('phi')
ax.set_zlabel('height')
ax.set_title('Planet')
fig.colorbar(surf, shrink=0.7, aspect=10) 
ax.view_init(90, 0)
             
plt.show()
