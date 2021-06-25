'''
Property of: NavINST Laboratory
Author: Marwan A. Rashed
Description: Test file of the mechanization class
'''
from datetime import time
from IntializationParameters import InitINS
from Mechanization import Mechanization
import pandas as pd 
import numpy as np
from scipy import io 
import math

import matplotlib.pyplot as plt

##### Read Data #####
###### IMU readings over time
IMU = pd.read_csv ('C:\\Users\\marwa\\Desktop\\INS\\INSfullMech\\Toronto_Traj_one_IMU_sync.csv')
# print (IMU)
###### Reference point 
Ref_Solution = pd.read_csv ('C:\\Users\\marwa\\Desktop\\INS\\INSfullMech\\Toronto_Traj_one_REF_sync.csv')

# print (Ref_Solution)

# print(Ref_Solution.iloc[0,3])
###### Initial position, attitudes
Init_lat, Init_long, Init_alt = Ref_Solution['Lat'][0], Ref_Solution['Long'][0], Ref_Solution['Alt'][0]
Init_roll, Init_pitch, Init_azimuth = Ref_Solution['Roll'][0], Ref_Solution['Pitch'][0], Ref_Solution['Azi'][0]

####### Take the specific imu readings needed as lists
timestamps =  list( IMU['Time Stamp'])
wx, wy , wz = list (IMU['w_x']),  list (IMU['w_y']), list (IMU['w_z'])
fx, fy, fz  = list (IMU['f_x']), list (IMU['f_y']), list (IMU['f_z'])
print(timestamps[0], wx[0], wy[0], wz[0], fx[0], fy[0], fz[0])

################################# Initialization Stage ########################
################# consruct an instance of the mechanization
INS_Mechanize = Mechanization (Init_lat, Init_long, Init_alt, Init_roll, Init_pitch, Init_azimuth)
INS_Mechanize.Init_Velocity(-0.0002, -0.0004, 0.0023)

################# Start iterating to get the INS solution ###############
delta_time = timestamps[1] - timestamps[0]
latitude, longitude, altitude = [0], [0], [0]
ve, vn , vu = [0], [0], [0]
roll, pitch, azimuth, yaw = [], [],[], []
g = 9.84
duration = len(timestamps) -1
latitude_prev,longitude_prev,altitude_prev = Init_lat, Init_long, Init_alt
ve_Prev, vn_Prev, vu_Prev = -0.0002, -0.0004, 0.0023

 
wx_bias = np.average(wx[0:1452])
wy_bias = np.average(wy[0:1452])
wz_bias = np.average(wz[0:1452])
fx_bias = np.average(fx[0:1452])
fy_bias = np.average(fy[0:1452])
fz_bias = np.average(fz[0:1452])

# wx -= wx_bias
# wy -= wy_bias
wz -= wz_bias
# fx -= fx_bias
# fy -= fy_bias
# fz -= fz_bias

print (wz_bias)
for i in range (1,5): 

    curr_latitude, curr_longitude, curr_altitude , curr_ve ,curr_vn, curr_vu, curr_roll, curr_pitch, curr_azimuth = INS_Mechanize.compile (0.2, 
                wx[i], wy[i], wz[i]
                , fx[i], fy[i],fz[i] 
                ,latitude_prev,longitude_prev,altitude_prev,
                 ve_Prev, vn_Prev, vu_Prev,
                 g ) 

    latitude.append(curr_latitude), longitude.append(curr_longitude), altitude.append(curr_altitude)
    ve.append(curr_ve), vn.append(curr_vn), vu.append(curr_vu)
    delta_time = timestamps [i] - timestamps[i-1]
    latitude_prev,longitude_prev,altitude_prev = curr_latitude, curr_longitude, curr_altitude
    ve_Prev, vn_Prev, vu_Prev = curr_ve ,curr_vn, curr_vu

plt.plot(latitude [1:-1], longitude [1:-1])
# plt.plot(Ref_Solution['Lat'].iloc[1:-1] , Ref_Solution['Long'].iloc[1:-1] )
plt.show()
