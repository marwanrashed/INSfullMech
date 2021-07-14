'''
Property of: NavINST Laboratory
Author: Marwan A. Rashed
Description: Validation file of the mechanization class on Kingston 2008 dataset
The validation will be run on: 
Frequency : 10 Hz 
Start time : 600 seconds
Biases :fx = -70.88e-4; fy = 40.80e-4; fz = 61.00e-4;      
        wx = -40.84e-007; wy = 12.3e-007; wz = 2.06e-006;

Note: fy and wy from novatel should be negative , hence an operation of:
    wy = - wy, and fy = -fy should always be proceeded in the beginning.
'''
import numpy as np
import pandas as pd
from Mechanization import Mechanization
from loading_dataset import syncronize_INS, synchronize_general, load_IMU_dataset, load_REF_dataset
import matplotlib.pyplot as plt
import HelperFuncs

####### files paths ########
ref_file_path = r'C:\Users\marwa\Desktop\INS\INSfullMech\2008-06-23 kingston downtown\rover\Mat\INSPVA'
imu_file_path = r'C:\Users\marwa\Desktop\INS\INSfullMech\2008-06-23 kingston downtown\rover\Mat\RAWIMU_denoised_LOD3_interpolated2'
Freq_INS = 10
KVH = 0
ref_data =  load_REF_dataset(ref_file_path, Freq_INS)
# Reference Keys are:  dict_keys(['time', 'lat', 'lon', 'alt', 'roll', 'pitch', 'azimuth', 'Ve', 'Vn', 'Vu'])
IMU_data = load_IMU_dataset(imu_file_path, KVH, Freq_INS)
# IMU Keys are:  dict_keys(['time', 'fx', 'fy', 'fz', 'wx', 'wy', 'wz'])
ref_data, IMU_data =  syncronize_INS(ref_data , IMU_data)

print ("The Imu original first timestamp", IMU_data['time'][0])
print ("The reference original first timestamp", ref_data['time'][0])
####### Declare biases #########
fx_bias = -70.88e-4 
fy_bias = 40.80e-4 
fz_bias = 61.00e-4      
wx_bias = -40.84e-007
wy_bias = 12.3e-007 
wz_bias = 2.06e-006
######## Remove Biases #########
IMU_data['fx'] -= fx_bias
IMU_data['fy'] -= fy_bias
IMU_data['fz'] -= fz_bias
IMU_data['wx'] -= wx_bias
IMU_data['wy'] -= wy_bias
IMU_data['wz'] -= wz_bias
############# Mechanization Intialization Params ###########
start_time = 600 # seconds
start = start_time * Freq_INS
Init_lat, Init_lon, Init_alt = ref_data ['lat'][start], ref_data ['lon'][start], ref_data ['alt'][start]
Init_roll, Init_pitch, Init_azimuth = ref_data ['roll'][start], ref_data ['pitch'][start], ref_data ['azimuth'][start]
Init_ve, Init_vn, Init_vu = ref_data ['Ve'][start], ref_data ['Vn'][start], ref_data ['Vu'][start]
dt = 0.1 # 10 Hz
############ Initialize Mecanization ##############
################# consruct an instance of the mechanization
INS_Mechanize = Mechanization (Init_lat, Init_lon, Init_alt, Init_roll, Init_pitch, Init_azimuth, dt)
INS_Mechanize.Init_Velocity(Init_ve, Init_vn, Init_vu)
duration = len(IMU_data['fx'])

Lat, Lon, Alt = [],[], []
Roll, Pitch, Azimuth = [],[],[]
ve,vn,vu = [],[],[]

################# Loop over the mechanization ##########
for i in range (start,duration): 
    wx, wy, wz = IMU_data['wx'][i], IMU_data['wy'][i], IMU_data['wz'][i]
    fx, fy, fz = IMU_data['fx'][i], IMU_data['fy'][i], IMU_data['fz'][i]
    INS_Mechanize.compile_standalone(wx, wy, wz,fx, fy, fz)
    Lat.append(INS_Mechanize.latitude)
    Lon.append(INS_Mechanize.longitude)
    Alt.append(INS_Mechanize.altitude)
    Roll.append(INS_Mechanize.roll)
    Pitch.append(INS_Mechanize.pitch)
    Azimuth.append(INS_Mechanize.azimuth)
    ve.append(INS_Mechanize.velocity_vector[0])
    vn.append(INS_Mechanize.velocity_vector[1])
    vu.append(INS_Mechanize.velocity_vector[2])

##############Prepare data for plotting##########
######## Prepare INS standalone and ref Output ######
imu_data_mechanized = dict()
ref_data_trimmed = dict()
imu_data_mechanized['time'], ref_data_trimmed['time'] = IMU_data['time'][start:duration], ref_data['time'][start:duration]
imu_data_mechanized['lat'], ref_data_trimmed['lat'] = Lat, ref_data['lat'][start:duration]
imu_data_mechanized['lon'], ref_data_trimmed['lon']  = Lon, ref_data['lon'][start:duration]
imu_data_mechanized['alt'], ref_data_trimmed['alt'] = Alt, ref_data['alt'][start:duration]
imu_data_mechanized['roll'], ref_data_trimmed['roll'] = Roll, ref_data['roll'][start:duration]
imu_data_mechanized['pitch'], ref_data_trimmed['pitch'] = Pitch, ref_data['pitch'][start:duration]
imu_data_mechanized['azimuth'], ref_data_trimmed['azimuth'] = Azimuth, ref_data['azimuth'][start:duration]
imu_data_mechanized['ve'], ref_data_trimmed['ve'] = ve, ref_data['Ve'][start:duration]
imu_data_mechanized['vn'], ref_data_trimmed['vn'] = vn, ref_data['Vn'][start:duration]
imu_data_mechanized['vu'], ref_data_trimmed['vu'] = vu, ref_data['Vu'][start:duration]

##############Plotting ###################
############## Full Diagnosis ################
print ("The Imu first timestamp", imu_data_mechanized['time'][0])
print ("The reference first timestamp", ref_data_trimmed['time'][0])
HelperFuncs.full_diagnosis ( ref_data_trimmed,  imu_data_mechanized, 'INS standalone')
