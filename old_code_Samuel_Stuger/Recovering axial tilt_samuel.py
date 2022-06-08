import matplotlib.pyplot as plt
import math
import numpy as np

'''
Parameters
'''
res = 15
rho = 6.371e6
R = 1.496e11
res_alpha = 37
res_beta = 19

alpha = np.linspace(0/4*np.pi,2/2*np.pi,res_alpha)
beta = np.linspace(0/4*np.pi,1/2*np.pi,res_beta)

day = 24
year = 365*day
omega_day = 2*np.pi/day
omega_year = 2*np.pi/year

delta_t = 1
delta_phi = np.pi/res
delta_theta = np.pi/res
delta_Omega = delta_phi*delta_theta

hours = np.arange(0,24)
time_array = hours.copy()

for i in range(1,20):
    time_array = np.append(time_array,hours+i*year/20)

time_res = len(time_array)

'''
Albedo-map
'''
A = np.zeros((2*res**2,res_alpha*res_beta))
f_arrays = np.zeros((480,res_alpha*res_beta))
time = np.linspace(0,1,round(year/delta_t))*365

'''
Eulerrotationmatrices
'''
def positive(arg):
    return(arg+abs(arg))/2

def y_rotation(angle):
    Y = np.array([[np.cos(angle),0,np.sin(angle)],[0,1,0],
                  [-np.s in (angle),0,np.cos(angle)]])
    return Y

def z_rotation(angle):
    Z = np.array([[np.cos(angle),-np.sin(angle),0],
                  [np.sin(angle),np.cos(angle),0],[0,0,1]])
    return Z

'''
Initialize transfromation matrices (edge-on, face-on)
'''

T = np.zeros((time_res,2*res**2,res_alpha*res_beta))
T_pinv = np.zeros((2*res**2,time_res,res_alpha*res_beta))

'''
Observer
'''
o_vec = np.array([0,0,1])

'''
Compute matrix elements
'''
for k in range(res_alpha):
    R_equinox = z_rotation(alpha[k])

    for l in range(res_beta):
        R_tilt = y_rotation(beta[l])
        R_axial = np.matmul(R_equinox,R_tilt)

        for i in range(len(time_array)):
            t = time_array[i]
            r_vec = np.array([-np.cos(omega_year*t),-np.sin(
            omega_year*t),0])
            R_daily = z_rotation(omega_day*t)
            daily_rotation = np.matmul(R_axial,R_daily)

            for j in range(2*res**2):
                phi = math.floor(j/res)*delta_phi
                theta = (j % res)*delta_theta

                s_vec = np.array([np.cos(phi)*np.sin(theta),np.sin(phi)*np.sin(theta),np.cos(theta)])
                s_vec_rotated = np.matmul(daily_rotation,s_vec)
                r_s = np.dot(r_vec,s_vec_rotated)
                s_o = np.dot(s_vec_rotated,o_vec)

                illuminated = positive(r_s)
                visible = positive(s_o)


                T[i][j][k*res_beta+l] = illuminated*visible*np.sin(theta)*delta_Omega

        T_ = T[:,:,k*res_beta+l]
        T_pinv_ = np.l in alg.pinv(T_,rcond = 0.01)
        T_pinv[:,:,k*res_beta+l] = T_pinv_
        A_ = np.matmul(T_pinv_,f_curve_edge)
        A[:,k*res_beta+l] = A_
        f_arrays[:,k*res_beta+l] = np.matmul(T_,A_)*rho**2/(R**2*np.pi)
        print(k*res_beta+l)

f_curve_matrix = np.tile(f_curve_edge,(res_alpha*res_beta,1)).transpose()
f_diff = np.square(f_curve_matrix-f_arrays)
square_diff = np.sum(f_diff,axis = 0)
minima = np.where(square_diff == square_diff.min())[0]

alpha_ = np.floor(minima/res_beta).astype( in t)
beta_ = minima % res_beta

print('alpha = ',np.degrees(alpha[alpha_]))
print('beta = ',np.degrees(beta[beta_]))

