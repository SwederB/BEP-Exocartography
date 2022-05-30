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
    
    theta_flat = np.ndarray.flatten(theta, order = 'C')
    phi_flat = np.ndarray.flatten(phi, order = 'C')
    
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
    

    coordinates = np.vstack((x,y,z)).transpose() # coordinates is (#points, 3) shaped

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


def polar(polarcaps, sphere,alts, p_ice):
    minz = min(sphere[:,2])
    maxz = max(sphere[:,2])
    rowdivide = np.linspace(minz,maxz,100)
    rangez = abs(maxz - minz)
    
    iceval = max(alts)

    for i in range(len(alts)):
        if abs(maxz-sphere[i,2])/rangez < p_ice or abs(minz-sphere[i,2])/rangez < p_ice:
            alts[i] = iceval
    return alts

'''
initialize sphere and final altitudes vector
'''
radius = 1
res = 50


sphere, theta, phi = create_sphere1(radius,res)
print('shape of sphere is', sphere.shape)
alts = np.zeros(len(sphere[:, 0]))
albedo_map = np.zeros(len(sphere[:, 0]))
albedo_map_comp = albedo_map.copy()
ocean_map = np.zeros(len(sphere[:, 0]))
ocean_map_comp = ocean_map.copy()
seeliger_map = np.zeros(len(sphere[:, 0]))
seeliger_map_comp = seeliger_map.copy()
'''
initialize vertices, seeds and altitudes(sea level)
'''


'''
calculate altitude for every point on sphere
'''
for i in range(len(sphere[:, 0])):
    vertices_ = np.array([[-2  ,  2.1,  2.2], [ 2  , -2.1,  2.2], [ 2  ,  2.1, -2.2],[-2.3,  -2.4,  -2.5]])
    seeds_ , altitudes_ = np.array([5,24,2,16]), np.array([0., 0. , 0., 0.])
    point = sphere[i, :]
    alt = height(point, vertices_, seeds_, altitudes_)
    alts[i] = alt
    print(i)
    
polarcaps = True
p_ice = 0.05

if polarcaps:
    alts = polar(polarcaps, sphere, alts, p_ice)
#    print(i, ': ', alt)

'''
albedo map values



'''


alt_grid1 = alts.reshape((res, 2*res), order = 'F')
alt_grid1_ = np.copy(alt_grid1)



alts1 = np.ndarray.flatten(alt_grid1_, order = 'F')




alsnow = 0.8
alsand = 0.4
alforest = 0.15
algrass = 0.25
alocean = 0.06

selsnow = 0.0
selsand = 0.5
selforest = 0.2
selgrass = 0.25
selocean = 0.01


'''
albedo map trial
'''

low20 = 0.2*max(alts)
high80 = 0.8*max(alts)


for i in range(len(sphere[:,0])):
    if alts1[i] < low20:
        albedo_map[i] = alocean
        ocean_map[i]  = 1
        seeliger_map[i] = selocean
    elif alts1[i] >= low20 and alts1[i] < 0.4*max(alts):
        albedo_map[i] = algrass
        seeliger_map[i] = selgrass
    elif alts1[i] >= 0.4*max(alts) and alts1[i] < 0.6*max(alts):
        albedo_map[i] = alforest
        seeliger_map[i] = selforest
    elif alts1[i] >= 0.6*max(alts) and alts1[i] < 0.8*max(alts):
        albedo_map[i] = alsand
        seeliger_map[i] = selsand
    else:
        albedo_map[i] = alsnow
        seeliger_map[i] = selsnow


print('busy plotting')
'''
create corresponding meshgrid
'''
alt_grid = alts1.reshape((res, 2*res), order = 'F')
albedo_map_grid = albedo_map.reshape((res, 2*res), order = 'F')
ocean_map_grid = ocean_map.reshape((res,2*res), order = 'F')
seeliger_map_grid = seeliger_map.reshape((res,2*res), order= 'F')

                                         
shrink = 0.5


newgistearth = cm.get_cmap('gist_earth', 4096)
newcolors = np.vstack((newgistearth(np.linspace(0,0.35,256)), newgistearth(np.linspace(0.45,1,256))))
newcmp = ListedColormap(newcolors)

albedo_colors = np.vstack((newgistearth(np.linspace(0.15,0.22,100)), newgistearth(np.linspace(0.45,0.55,200)),newgistearth(np.linspace(0.6,0.7,200)),newgistearth(np.linspace(0.7,0.8,200)), newgistearth(np.linspace(0.95,1,400))))
albedo_cmp = ListedColormap(albedo_colors)
albedo_colors1 = np.vstack((newgistearth(np.linspace(0.15,0.22,140)), newgistearth(np.linspace(0.45,0.55,180)),newgistearth(np.linspace(0.6,0.7,200)),newgistearth(np.linspace(0.7,0.8,280)), newgistearth(np.linspace(0.95,1,400))))
albedo_cmp1 = ListedColormap(albedo_colors1)

lon = np.linspace(-np.pi,np.pi,2*res)
lat = np.linspace(-np.pi/2,np.pi/2,res)
Lon,Lat = np.meshgrid(lon,lat)

fig1 = plt.figure(1, figsize = (20,15))
ax1 = fig1.add_subplot(131, projection = 'mollweide')
ax2 = fig1.add_subplot(132, projection = 'mollweide')
ax3 = fig1.add_subplot(133, projection = 'mollweide')

ax1.title.set_text('Lambertian Albedo')
ax2.title.set_text('Ocean "Albedo"')
ax3.title.set_text('Seeliger "Albedo"')


im1 = ax1.pcolormesh(Lon,Lat, np.flipud(albedo_map_grid),cmap = albedo_cmp1)
plt.colorbar(im1,shrink = shrink, aspect=40,ax = ax1)
tick_labels = np.array(['-150°', '', '-90°','', '-30°', '', '30°', '','90°', '', '150°'])
tick_labels_y = ['-75°', '', '-45°','', '-15°', '', '15°', '','45°', '', '75°']


im2 = ax2.pcolormesh(Lon,Lat, np.flipud(ocean_map_grid),cmap = 'binary')
plt.colorbar(im2,shrink = shrink, aspect=40,ax = ax2)
tick_labels = np.array(['-150°', '', '-90°','', '-30°', '', '30°', '','90°', '', '150°'])
tick_labels_y = ['-75°', '', '-45°','', '-15°', '', '15°', '','45°', '', '75°']


     
im3 = ax3.pcolormesh(Lon,Lat, np.flipud(seeliger_map_grid),cmap = 'Greys')
plt.colorbar(im3,shrink = shrink, aspect=40, ax = ax3)
tick_labels = np.array(['-150°', '', '-90°','', '-30°', '', '30°', '','90°', '', '150°'])
tick_labels_y = ['-75°', '', '-45°','', '-15°', '', '15°', '','45°', '', '75°']
                                
                                         
             
plt.show()
