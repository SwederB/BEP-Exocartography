import matplotlib.pyplot as plt
import math
import numpy as np

'''
Parameters
'''
res = 45
rho = 6.371e6
R = 1.496e11
alpha = 60/180*np.pi
beta = 60/180*np.pi

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

N_ave = 187

for i in range(1,200):
    time_array = np.append(time_array,hours+i*year/200)

time_res = len(time_array)

'''
Albedo-map
'''
A = np.ones(2*res**2)
time = np.linspace(0,1,round(year/delta_t))*365


'''
Euler rotation matrices
'''
def positive(arg):
    return(arg+abs(arg))/2

def y_rotation(angle):
    Y = np.array([[np.cos(angle),0,np.sin(angle)],[0,1,0],
                  [-np.sin(angle),0,np.cos(angle)]])
    return Y

def z_rotation(angle):
    Z = np.array([[np.cos(angle),-np.sin(angle),0],
                  [np.sin(angle),np.cos(angle),0],[0,0,1]])
    return Z



R_equinox = z_rotation(alpha)
R_tilt = y_rotation(beta)
R_axial = np.matmul(R_equinox,R_tilt)


'''
Initialize transfromation matrices (edge-on, face-on)
'''
T = np.zeros((time_res,2*res**2))

T1 = T.copy()

'''
Observer
'''
o_vec = np.array([1,0,0])
o_vec_1 = np.array([0,0,1])


'''
Compute matrix elements
'''
for i in range(len(time_array)):
    t = time_array[i]
    r_vec = np.array([-np.cos(omega_year*t),-np.sin(omega_year*t),0])
    R_daily = z_rotation(omega_day*t)
    daily_rotation = np.matmul(R_axial,R_daily)

    for j in range(2*res**2):
        phi = math.floor(j/res)*delta_phi
        theta = (j % res)*delta_theta

        s_vec = np.array([np.cos(phi)*np.sin(theta),
                          np.sin(phi)*np.sin(theta),np.cos(theta)])
        s_vec_rotated = np.matmul(daily_rotation,s_vec)

        r_s = np.dot(r_vec,s_vec_rotated)
        s_o = np.dot(s_vec_rotated,o_vec)
        s_o_1 = np.dot(s_vec_rotated,o_vec_1)


        illuminated = positive(r_s)
        visible = positive(s_o)
        visible_1 = positive(s_o_1)

        T[i][j] = illuminated*visible*np.sin(theta)*delta_Omega
        T1[i][j] = illuminated*visible_1*np.sin(theta)*delta_Omega



'''
Compute light-curves
'''
f_curve_edge = np.matmul(T,albedo_map)*rho**2/(R**2*np.pi)
f_curve_face = np.matmul(T1,albedo_map)*rho**2/(R**2*np.pi)


noise_poisson = np.random.poisson(N_ave,time_res)/N_ave
noise_gauss = np.random.normal(0,1,time_res)

#f_curve_edge_noisy = f_curve_edge+noise_gauss*max(f_curve_edge)/100
#f_curve_face_noisy = f_curve_face+noise_gauss*max(f_curve_face)/100

f_curve_edge_noisy = np.multiply(f_curve_edge,noise_poisson)
f_curve_face_noisy = np.multiply(f_curve_face,noise_poisson)

'''
Save transfer matrices for later use
'''
#np.savez_compressed('Albedo_a90_b90_20d_res_15.npz',edge = T,face = T1)


'''
Compute and save SVD
'''
SVD_edge = np.l in alg.svd(T,compute_uv = False)
SVD_face = np.l in alg.svd(T1,compute_uv = False)
'''
Compute inverse transfer matrix
'''
T_pinv_edge = np.linalg.pinv(T,rcond = 0.025)
T_pinv_face = np.linalg.pinv(T1,rcond = 0.025)
A_edge = np.matmul(T_pinv_edge,f_curve_edge)
A_face = np.matmul(T_pinv_face,f_curve_face)


edge_factor = 0.8/(max(A_edge)-min(A_edge))
face_factor = 0.8/(max(A_face)-min(A_face))


A_edge_scaled = (np.matmul(T_p in v_edge,f_curve_edge)-m in (A_edge))*edge_factor
A_face_scaled = (np.matmul(T_p in v_face,f_curve_face)-m in (A_face))*face_factor

reconstuct_A_grid_edge_scaled = A_edge_scaled.reshape((res,2*res),order = 'F')
reconstuct_A_grid_face_scaled = A_face_scaled.reshape((res,2*res),order = 'F')

edge_diff_grid = np.abs(albedo_map_grid-reconstuct_A_grid_edge_scaled)
face_diff_grid = np.abs(albedo_map_grid-reconstuct_A_grid_face_scaled)

theta_array_plus = np.linspace(0,np.pi,res)+1/2*delta_theta
theta_array_min = np.linspace(0,np.pi,res)-1/2*delta_theta

theta_array_plus[res-1] = np.pi
theta_array_min[0] = 0

theta_matrix_plus = np.tile(theta_array_plus,(2*res,1)).transpose()
theta_matrix_min = np.tile(theta_array_min,(2*res,1)).transpose()

surf_matrix = delta_phi*(np.cos(theta_matrix_min)-np.cos(theta_matrix_plus))/(4*np.pi)

edge_error = np.sum(np.multiply(np.square(edge_diff_grid),surf_matrix))
face_error = np.sum(np.multiply(np.square(face_diff_grid),surf_matrix))

# print(edge_error)
# print(face_error)



A_edge_noisy = np.matmul(T_pinv_edge, f_curve_edge_noisy)
A_face_noisy = np.matmul(T_pinv_face, f_curve_face_noisy)


edge_factor_noisy = 0.8/(max(A_edge_noisy)-min(A_edge_noisy))
face_factor_noisy = 0.8/(max(A_face_noisy)-min(A_face_noisy))


A_edge_scaled_noisy = A_edge_noisy*edge_factor_noisy
A_face_scaled_noisy = A_face_noisy*face_factor_noisy

reconstuct_A_grid_edge_scaled_noisy = A_edge_scaled_noisy.reshape((res,2*res),order = 'F')
reconstuct_A_grid_face_scaled_noisy = A_face_scaled_noisy.reshape((res,2*res),order = 'F')

edge_diff_grid_noisy = np.abs(albedo_map_grid-reconstuct_A_grid_edge_scaled_noisy)
face_diff_grid_noisy = np.abs(albedo_map_grid-reconstuct_A_grid_face_scaled_noisy)


edge_error_noisy = np.sum(np.multiply(np.square(edge_diff_grid_noisy),surf_matrix))
face_error_noisy = np.sum(np.multiply(np.square(face_diff_grid_noisy),surf_matrix))

print(edge_error_noisy)
print(face_error_noisy)

