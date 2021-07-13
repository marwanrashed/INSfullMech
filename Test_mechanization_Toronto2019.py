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
start_time =  1600
###### Initial position, attitudes
Init_lat, Init_long, Init_alt = Ref_Solution['Lat'][start_time], Ref_Solution['Long'][start_time], Ref_Solution['Alt'][start_time]
Init_roll, Init_pitch, Init_azimuth = Ref_Solution['Roll'][start_time], Ref_Solution['Pitch'][start_time], Ref_Solution['Azi'][start_time]

####### Take the specific imu readings needed as lists
timestamps =  list( IMU['Time Stamp'])
wx, wy , wz = list (IMU['w_x']),  list (IMU['w_y']), list (IMU['w_z'])
fx, fy, fz  = list (IMU['f_x']), list (IMU['f_y']), list (IMU['f_z'])
Ve, Vn, Vu = Ref_Solution['Ve'], Ref_Solution['Vn'], Ref_Solution['Vu']
ROLL, PITCH, AZIMUTH =  Ref_Solution['Roll'], Ref_Solution['Pitch'], Ref_Solution['Azi']
LAT, LON, ALT = Ref_Solution['Lat'], Ref_Solution['Long'], Ref_Solution['Alt']
print(timestamps[0], wx[0], wy[0], wz[0], fx[0], fy[0], fz[0])
print (Init_lat, Init_long, Init_alt, Init_roll, Init_pitch, Init_azimuth)
################################# Initialization Stage ########################
################# consruct an instance of the mechanization
INS_Mechanize = Mechanization (Init_lat, Init_long, Init_alt, Init_roll, Init_pitch, Init_azimuth, 0.2)
INS_Mechanize.Init_Velocity(Ve[start_time], Vn[start_time], Vu[start_time])

################# Start iterating to get the INS solution ###############
delta_time = timestamps[start_time+1] - timestamps[start_time]
latitude, longitude, altitude = [Init_lat], [Init_long], [Init_alt]
ve, vn , vu = [-0.0002], [-0.0004], [0.0023]
roll, pitch, azimuth = [Init_roll], [Init_pitch],[Init_azimuth]
g = 9.84
duration = len(timestamps) -1
latitude_prev,longitude_prev,altitude_prev = Init_lat, Init_long, Init_alt
ve_Prev, vn_Prev, vu_Prev = Ve[start_time], Vn[start_time], Vu[start_time]
roll_prev, pitch_prev, azimuth_prev = Init_roll, Init_pitch, Init_azimuth

print (duration)
bias_time = 50
wx_bias = np.mean(wx[0:bias_time])
wy_bias = np.mean(wy[0:bias_time])
wz_bias =  - np.mean(wz[0:bias_time]) #-0.001565
fx_bias =   np.mean(fx[0:bias_time]) #0.0725
fy_bias =  np.mean(fy[0:bias_time]) #0.2258
fz_bias = np.mean(fz[0:bias_time])

wx -= wx_bias
wy -= wy_bias
wz -= wz_bias #list(map(lambda x : x- wz_bias , wz )) 
# fx = #list(map(lambda x : x- fx_bias  , fx )) 
# fy = #list (map(lambda x : x- fy_bias , fy )) 
fx -= fx_bias
fy -= fy_bias 
fz -= fz_bias
duration_2 = start_time + 300
# print (fx_bias, fy_bias)
for i in range (start_time+1 , duration_2 ): 
    # print ("############################################### Iteration NO % ###############################", i)
    # print ("Reference Attitudes", ROLL[i], PITCH[i] )
    # print ("Reference azimuth", AZIMUTH[i])
    # print ("Reference Velocities", [Ve[i], Vn[i], Vu[i]])
    # print ("Reference Position",[LAT[i], LON[i], ALT[i]])
    # INS_Mechanize.compile_standalone ( 
    #             wx[i], wy[i], wz[i]
    #             , fx[i], fy[i],fz[i]) 
    
    # INS_Mechanize.compile_closed_loop ( 
    #         wx[i], wy[i], wz[i]
    #         , fx[i], fy[i],fz[i]
    #         , latitude_prev,longitude_prev,altitude_prev
    #         , [ve_Prev, vn_Prev, vu_Prev]
    #         , roll_prev, pitch_prev, azimuth_prev ) 
    curr_latitude, curr_longitude, curr_altitude = INS_Mechanize.latitude , INS_Mechanize.longitude, INS_Mechanize.altitude
    curr_ve ,curr_vn, curr_vu = INS_Mechanize.velocity_vector
    curr_roll, curr_pitch, curr_azimuth = INS_Mechanize.roll, INS_Mechanize.pitch, INS_Mechanize.azimuth
    latitude.append(curr_latitude), longitude.append(curr_longitude), altitude.append(curr_altitude)
    ve.append(curr_ve), vn.append(curr_vn), vu.append(curr_vu)
    delta_time = timestamps [i] - timestamps[i-1]
    latitude_prev,longitude_prev,altitude_prev = curr_latitude, curr_longitude, curr_altitude
    ve_Prev, vn_Prev, vu_Prev = curr_ve ,curr_vn, curr_vu
    roll_prev, pitch_prev, azimuth_prev = curr_roll, curr_pitch, curr_azimuth

plt.plot(longitude  [0:-1], latitude [0:-1])
plt.plot(Ref_Solution['Long'].iloc[start_time:duration_2], Ref_Solution['Lat'].iloc[start_time:duration_2])
plt.show()
