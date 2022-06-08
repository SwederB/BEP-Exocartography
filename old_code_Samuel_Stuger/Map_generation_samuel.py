import matplotlib.pyplot as plt
import numpy as np
import random
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay , ConvexHull
from matplotlib import cm
from matplotlib.colors import ListedColormap ,LinearSegmentedColormap

def create_sphere1(r,res): #r = radius, res = resolution
    '''
    creates points on a sphere in the form :
    [ [ x1 . y1 . z1 . ]
    . . . .
    [ xn . yn . zn . ] ]
    '''
    phi = np.linspace(0, 2*np.pi, 2*res)
    theta = np.linspace(0, np.pi, res)
    
    theta , phi = np.meshgrid(theta, phi)
    theta_flat = np.ndarray.flatten(theta, order = 'C')
    phi_flat = np.ndarray.flatten(phi, order = 'C')

    r_pre = r*np.sin(theta_flat)
    x = np.cos(phi_flat)*r_pre
    y = np.sin(phi_flat)*r_pre
    z = r*np.cos(theta_flat)

    coordinates = np.vstack((x, y, z)).transpose()

    return coordinates, theta , phi




def pnt_in_cvex_hull(hull,pnt): 
    '''checks if 'pnt' is inside the convex hull'''
    new_hull = ConvexHull(np.concatenate((hull.points, [pnt])))
    if np.array_equal(new_hull.vertices, hull.vertices):
        return True
    else:
        return False 

    
    
def height(p,vertices,seeds,altitudes):
    max_dist=0
    '''
    finding the shortest edge
    '''
    for i in range(3):
        for j in range(i+1,4):
            dist = np.linalg.norm(vertices[i]-vertices[j])
            if dist > max_dist:
                max_dist = dist
                furthest = np.array([i,j])

    '''
    reordering everything
    '''
    vertices[[0,furthest[0]]] = vertices[[furthest[0],0]]
    vertices[[1,furthest[1]]] = vertices[[furthest[1],1]]

    seeds[[0,furthest[0]]] = seeds[[furthest[0],0]]
    seeds[[1,furthest[1]]] = seeds[[furthest[1],1]]

    altitudes[[0,furthest[0]]] = altitudes[[furthest[0],0]]
    altitudes[[1,furthest[1]]] = altitudes[[furthest[1],1]]

    '''
    creating new edge
    newaltitude = average altitude + random(-0,05;0,05)*sqrt(
    distance)
    '''
    v_new = (vertices[0]+vertices[1])/2
    s_new = (seeds[0]+seeds[1])/2
    random.seed(s_new)
    a_new = (altitudes[0]+altitudes[1])/2+0.01*(random.random()
    -0.5)*max_dist**1.5

    '''
    finding in which tetrahedon the point p is
    '''
    tetra = np.copy(vertices)
    tetra[1] = v_new
    hull = ConvexHull(tetra)

    if pnt_in_cvex_hull(hull,p):
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
    if max_dist <= 0.001:
        return np.sum(altitudes)/4
    else:
        height(p,vertices,seeds,altitudes)
    return np.sum(altitudes)/4




'''
initialize sphere and final altitudes vector
'''
res = 45
sphere,theta,phi = create_sphere1(1,res)
alts = np.zeros(len(sphere[:,0]))
albedo_map = np.zeros(len(sphere[:,0]))

albedo_map_comp = albedo_map.copy()

'''
initialize vertices, seeds and altitudes(sea level)



calculate altitude for every point on sphere
'''
for i in range(len(sphere[:,0])):
    vertices_ = np.array([[-2,2.1,2.2],[2,-2.1,2.2],[2,2.1,-2.2],[-2.3,-2.4,-2.5]])

    seeds_,altitudes_ = np.array([1,2,3,4]),np.array([0.,0.,0.,0.])
    point = sphere[i,:]
    alt = height(point,vertices_,seeds_,altitudes_)
    alts[i] = alt

    
'''
calculate albedomap
'''

ave = np.sum(alts)/len(alts)
alts = alts-ave
max_alt = max(alts)
min_alt = min(alts)

alt_grid1 = alts.reshape((res,2*res),order = 'F')
alt_grid1_ = np.copy(alt_grid1)

alts1 = np.ndarray.flatten(alt_grid1_,order = 'F')


for i in range(len(sphere[:,0])):
    if alts1[i] >= 0.8*max_alt:
        albedo_map[i] = 0.8 #snow
    elif alts1[i] < 0.0:
        albedo_map[i] = 0.06 #ocean
    elif (alts1[i] >= 0.4*max_alt and alts1[i]<0.8*max_alt):
        albedo_map[i] = 0.4 #soil
    elif (alts1[i] >= 0.0 and alts1[i]<0.4*max_alt):
        albedo_map[i] = 0.15 #forest

'''
create corresponding meshgrid
'''
alt_grid = alts1.reshape((res,2*res),order = 'F')
albedo_map_grid = albedo_map.reshape((res,2*res),order = 'F')

new_gist_earth = cm.get_cmap('gist_earth',4096)
newcolors = np.vstack((new_gist_earth(np.linspace(0,0.35,256)), 
                       new_gist_earth(np.linspace(0.45,1,256))))
newcmp = ListedColormap(newcolors)
albedo_colors = np.vstack((new_gist_earth(np.linspace(0.15,0.22,100)),
                           new_gist_earth(np.linspace(0.45,0.55,200)),
                           new_gist_earth(np.linspace(0.7,0.8,300)),
                           new_gist_earth(np.linspace(0.95,1,400))))
albedo_cmp = ListedColormap(albedo_colors)
albedo_colors1 = np.vstack((new_gist_earth(np.linspace(0.15,0.22,140)),
                            new_gist_earth(np.linspace(0.45,0.55,180)),
                            new_gist_earth(np.linspace(0.7,0.8,280)),
                            new_gist_earth(np.linspace(0.95,1,400))))





albedo_cmp1 = ListedColormap(albedo_colors1)
lon = np.linspace(-np.pi,np.pi,2*res)
lat = np.linspace(-np.pi/2.,np.pi/2,res)
Lon,Lat = np.meshgrid(lon,lat)

plotalts= np.flipud(alt_grid)
# try to rotate the plot by 90 degrees by moving the central meridian 90 degrees
plotalts = np.roll(plotalts,plotalts.shape[1]//4*3,axis=1)

fig1 = plt.figure(1,figsize = (20,9))
ax1 = fig1.add_subplot(111)
im1 = ax1.pcolormesh(Lon,Lat,plotalts,cmap = newcmp)
plt.colorbar(im1,shrink = 0.75,aspect = 40)
# fig1.show()

