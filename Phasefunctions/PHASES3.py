import matplotlib.ticker as tck
import matplotlib.pyplot as plt
import numpy as np
import random
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay, ConvexHull
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import math



def observabledomain(cos_theta_i,cos_theta_r):
    if cos_theta_i >= 0 and cos_theta_r >= 0:
        return 1
    else:
        return 0    
def deltapeak(x,sigma):
    delta = np.exp(-(x)**2/(2*sigma))/np.sqrt(2*np.pi*sigma)
    return delta

def format_func(value, tick_number):
    # find number of multiples of pi/6
    N = int(np.round(6 * value / np.pi))
    if N == 0:
        return "0"
    elif N == 1:
        return r"$\pi/6$"
    elif N == 2:
        return r"$2\pi/6$"
    elif N == 3:
        return r"$3\pi/6$"
    elif N == 4:
        return r"$4\pi/6$"
    elif N == 5:
        return r"$5\pi/6$"
    elif N == 6:
        return r"$\pi$"
    elif N % 6 > 0:
        return r"${0}\pi/6$".format(N)
    else:
        return r"${0}\pi$".format(N // 6)



res = 20
N_phi = res
N_theta = res
N_alpha = 100


sigma = 0.01
nwater = 1.33
epsilon = 0.000000001
Lamb = np.zeros((N_alpha,2*res**2))
Lambcheck = np.zeros((N_alpha,2*res**2))
Opp1 = np.zeros((N_alpha,2*res**2))
Opp2 = np.zeros((N_alpha,2*res**2))
Lommel = np.zeros((N_alpha,2*res**2))
normcheck = np.zeros((N_alpha,2*res**2))




delta_phi = np.pi/res
delta_theta = np.pi/res
delta_Omega = delta_phi*delta_theta




alpha = np.linspace(0,np.pi,N_alpha)
m = 0

for i in range(N_alpha):
    i_hat = np.array([np.cos(alpha[i]),np.sin(alpha[i]),0])
    r_hat = np.array([1,0,0])
    print(str(i+1)+' /'+str(N_alpha))

    for j in range(2*res**2):
        phi = math.floor(j/res)*delta_phi
        theta = (j%res)*delta_theta
        s_hat = np.array([np.cos(phi)*np.sin(theta),np.sin(phi)*np.sin(theta),np.cos(theta)])


        cos_theta_i = np.inner(i_hat,s_hat)
        theta_i = np.arccos(cos_theta_i)
        

        cos_theta_r = np.inner(r_hat,s_hat)
        theta_r = np.arccos(cos_theta_r)


##        p1 = np.cross(i_hat,s_hat)/(np.linalg.norm(np.cross(i_hat,s_hat))+epsilon) #epsilon to ensure norm(p1)<1
##        
##
##        p2 = np.cross(r_hat,s_hat)/(np.linalg.norm(np.cross(r_hat,s_hat))+epsilon)
##        
##
##
##        cos_psi = np.inner(p1,p2)
##        
##        
##        psi = np.arccos(cos_psi)
##
##        Rs = abs((np.cos(theta_i)-np.sqrt(nwater**2-np.sin(theta_r)**2))/(np.cos(theta_i)+np.sqrt(nwater**2-np.sin(theta_i)**2)))**2
##        Rp = abs(((nwater**2)*np.cos(theta_i)-np.sqrt(nwater**2-np.sin(theta_r)**2))/((nwater**2)*np.cos(theta_i)+np.sqrt(nwater**2-np.sin(theta_i)**2)))**2
##
        Lamb[i][j]    = 4*observabledomain(cos_theta_i,cos_theta_r)*cos_theta_i*cos_theta_r*np.sin(theta)*delta_Omega/np.pi
##        normcheck[i][j] = np.sin(theta)*delta_Omega
##        Opp1[i][j]    = 4*observabledomain(cos_theta_i,cos_theta_r)*cos_theta_i*cos_theta_r*np.sin(theta)*delta_Omega*deltapeak(theta_i-theta_r,sigma)*deltapeak(psi-np.pi,sigma)/(np.sin(theta_i)*np.cos(theta_i)+epsilon)
##        Opp2[i][j]    = Opp1[i][j]*(Rs+Rp)

        Lommel[i][j]  = 2/(3*np.pi**2)*4*observabledomain(cos_theta_i,cos_theta_r)*cos_theta_i*cos_theta_r*np.sin(theta)*delta_Omega*(np.sin(alpha[i])+(np.pi-alpha[i])*np.cos(alpha[i]))/(np.cos(theta_i)+np.cos(theta_r)+epsilon)
        

        
A = np.ones(2*res**2)

lamb_curve = np.matmul(Lamb,A)
normcheck_curve = np.matmul(normcheck,A)
lambertian_analytical = 8/(3*np.pi)*(np.sin(alpha)+(np.pi-alpha)*np.cos(alpha))
##
##
##opp1_curve = np.matmul(Opp1,A)
##opp2_curve = np.matmul(Opp2,A)
##
##wp = 3
##
##opp1_curve[:wp] = [1]*wp
##opp2_curve[:wp] = [0]*wp
lommel_curve = np.matmul(Lommel,A)
##
##


f,ax=plt.subplots(figsize=(10,5))
ax.plot(alpha,normcheck_curve, label= 'Unity')
ax.plot(alpha,lamb_curve, label= 'Lambertian')
ax.plot(alpha, lambertian_analytical, label= 'Analytical')
ax.plot(alpha,lamb_curve/2, label= 'Lambertian/2')


##ax.plot(alpha, opp1_curve, label= 'Metallic glint')
##ax.plot(alpha, opp2_curve, label= 'Water glint (Fresnel)')
##ax.plot(alpha, lommel_curve, label= 'Lommel-Seeliger')
plt.xlabel('alpha')
plt.ylabel('f_curve')
plt.title('Phase function')
plt.legend()
plt.grid()
ax.xaxis.set_major_locator(plt.MultipleLocator(np.pi / 6))
ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 12))
ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))


plt.show()
