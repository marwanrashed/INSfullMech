import numpy as np
import pandas as pd
from Mechanization import Mechanization
from loading_dataset import syncronize_INS, synchronize_general, load_IMU_dataset, load_REF_dataset, load_IMU_data_no_downsampling
import matplotlib.pyplot as plt
import HelperFuncs
import scipy.io as sio



####### files paths ########
# ref_file_path = r'C:\Users\marwa\Desktop\INS\INSfullMech\2008-06-23 kingston downtown\rover\Mat\INSPVA'
imu_or_file_path = r'C:\Users\marwa\Desktop\INS\INSfullMech\2008-06-23 kingston downtown\rover\Mat\RAWIMU_denoised_LOD3_interpolated2'

imu_file_path = r'C:\Users\marwa\Desktop\INS\INSfullMech\f_w_IMSseconds_dt_Marwan.mat'
sample_file_path = r'C:\Users\marwa\Desktop\INS\INSfullMech\sample_time.mat'

sample_time_new = sio.loadmat(sample_file_path)
sample_time_new_list = sample_time_new['sample_time'][0][:]
print (len(sample_time_new_list))
print(min (sample_time_new_list))
print(max (sample_time_new_list))

imu_data_or = load_IMU_dataset (imu_or_file_path, 0, 10)
imu_data = load_IMU_data_no_downsampling (imu_file_path,0)
print ('Keys of IMU data', imu_data.keys())

print ('IMU Data:', imu_data)
print ('IMU Data fy sign:', imu_data['wy'][6000], 'IMU data original ',imu_data_or['wy'][6000] )
print (len(imu_data['wy']))
# imu_data = sio.loadmat(imu_file_path)
# #get all important data
# #Time
# time_imu = imu_data.get('IMU_second')[0][:] 
# print (time_imu)
# #sample time
# dt_imu = imu_data.get('sample_time')[0][:]
# print (dt_imu)

# #get imu frequency
# Freq_imu = np.round(1 / np.mean(time_imu[1:] - time_imu[:-1]))
# print('imu frequency = ' , Freq_imu , ' Hz')
# #-----------------------------------------------
# #accelorometer data
# fx_imu = imu_data.get('f')['x'][0,0][0,:]
# print (fx_imu)

# fy_imu = imu_data['f']['y'][0,0][0,:]
# print (fy_imu)
# # if IMU_CHOICE == KVH:
# #     fy_imu = fy_imu* -1 # special case for the accelerometer of imu is inverted
# fz_imu = imu_data['f']['z'][0,0][0,:]
# print (fz_imu)

# #-----------------------------------------------
# #Gyro data
# wx_imu = imu_data['w']['x'][0,0][0,:]
# wy_imu = imu_data['w']['y'][0,0][0,:] 
# # if IMU_CHOICE == KVH:
# #     wy_imu = wy_imu* -1 # special case for the accelerometer of imu is inverted
# wz_imu = imu_data['w']['z'][0,0][0,:] 


# imu_data = sio.loadmat(imu_file_path)

# ref_data =  load_REF_dataset(ref_file_path, Freq_INS)
# # Reference Keys are:  dict_keys(['time', 'lat', 'lon', 'alt', 'roll', 'pitch', 'azimuth', 'Ve', 'Vn', 'Vu'])
# IMU_data = load_IMU_dataset(imu_file_path, KVH, Freq_INS)
# # IMU Keys are:  dict_keys(['time', 'fx', 'fy', 'fz', 'wx', 'wy', 'wz'])
# ref_data, IMU_data =  syncronize_INS(ref_data , IMU_data)

# print ("The Imu original first timestamp", IMU_data['time'][0])
# print ("The reference original first timestamp", ref_data['time'][0])

# IMU_dataframe = pd.DataFrame (IMU_data)

# print ("IMU dataframe", IMU_dataframe.iloc[6000:]['time'].head(10))