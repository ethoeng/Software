import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from math import floor

def UpdateVelocityBoris(q,m,v_vec,B_vec,dt_sec):
    # print('v_in:',v_vec,'\t')
    t_vec = np.multiply((q/m)*(dt_sec/2),B_vec)
    s_vec = np.multiply(2/(1+np.linalg.norm(t_vec)**2),t_vec)
    # print('t_vec',t_vec, '\t','s_vec:',s_vec,'\t')

    v_minus_vec = v_vec #+ (q/m)*(dt_sec/2)*E_vec
   
    v_prime_vec = v_minus_vec + np.cross(v_minus_vec,t_vec)
    # print('v_prime:',v_prime_vec,'\t')

    v_plus_vec = v_minus_vec + np.cross(v_prime_vec,s_vec)
    # print('v_out:',v_plus_vec,'\n')

    return v_plus_vec #+ (q/m)*(dt_sec/2)*E_vec


def PushParticle(r_vec, v_vec, dt_sec):
    r_vec[0] += np.multiply(v_vec[0],dt_sec)
    r_vec[1] += np.multiply(v_vec[1],dt_sec)
    r_vec[2] += np.multiply(v_vec[2],dt_sec)
    return r_vec

"""
Initiate Fields and Particles
"""
    
# Load Fields (and interpolate)
filename = '2kG_CST.txt'#'2kG_Bn.txt'
Bx_T_data = np.loadtxt(filename,skiprows=2,usecols=1)#np.loadtxt(filename,comments='%',usecols=1)
r_m_data = np.loadtxt(filename,skiprows=2,usecols=0)*1e-3#np.loadtxt(filename,comments='%',usecols=0)*1e-3
Bx_interp_func = interpolate.interp1d(r_m_data,Bx_T_data,fill_value='extrapolate')

# Initiate Particle (SampleParticle)
ion_mass_amu = 7
amu_to_kg = 1.67e-27
m_kg = ion_mass_amu*amu_to_kg

ion_charge = +1
ion_charge_C = ion_charge*1.6e-19

KE_init_eV = 20000
eV_to_joule = 1.6e-19
KE_init_joule = KE_init_eV*eV_to_joule


y0 = 0
z0 = 1.5
vy_0 = 0
vz_0 = -np.sqrt(2*KE_init_joule/(m_kg))

r_vec=[0,y0,z0]
v_vec=[0,vy_0,vz_0]
print('r_init:',r_vec,' v_init: ',v_vec)

y_arr = [y0]
z_arr = [z0]
vy_arr=[vy_0]
vz_arr=[vz_0]

# Grid size
dr_m = np.abs(r_m_data[1]-r_m_data[0])
dt_sec = np.abs(dr_m/vz_0)
it_max = floor(r_m_data[-1]/dr_m)
print(dt_sec)

""" Push velocity back in time by 1/2 dt """
# E = EvalE(part.x)
# B = EvalB(part.x)

B_vec = [Bx_interp_func(np.sqrt(r_vec[1]**2+r_vec[2]**2)),0,0]
# B_vec = [0,0,0.01]
B_arr = [B_vec[0]]

# UpdateVelocity(part,E,B,-0.5*dt)
v_vec = UpdateVelocityBoris(ion_charge_C,m_kg,v_vec,B_vec,-0.5*dt_sec)

print('v_pushed_1/2',v_vec)

for it in range(it_max):
#B=EvalB(part.x)
    # print(it,r_vec)
    B_vec = [Bx_interp_func(np.sqrt(r_vec[1]**2+r_vec[2]**2)),0,0]
    # B_vec = [1,0,0]
#UpdateVelocity(part,E,B,dt)
    v_vec = UpdateVelocityBoris(ion_charge_C,m_kg,v_vec, B_vec,dt_sec)
    

#PushParticle(part,dt)
    r_vec = PushParticle(r_vec,v_vec,dt_sec)
    
#Store Data for plot
    B_arr.append(B_vec[0])
    
    vy_arr.append(v_vec[1])
    vz_arr.append(v_vec[2])

    y_arr.append(r_vec[1])
    z_arr.append(r_vec[2])
    
    # print(B_vec, '\t',v_vec,'\t',r_vec)

plt.plot(z_arr,np.multiply(y_arr,1e3),'*')
plt.xlabel('z [m]')
plt.ylabel('deflection [mm]')
plt.ylim(bottom=-25,top=30)
plt.yticks(np.arange(-25, 30, step=5.0))
plt.grid('on')

plt.plot(z_arr,B_arr,'*')