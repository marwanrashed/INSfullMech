import os
import scipy.io as sio
import numpy as np
#################################################################################
######################### Helper Functions ######################################
#################################################################################
def down_sample(signal , down_sample_bin ):
    L = len(signal)
    L_new = int(np.floor(L/down_sample_bin))
    signal_new = np.zeros(L_new)

    for i in range(L_new):
        frm = i * down_sample_bin
        to = (i+1) * down_sample_bin
        signal_new[i] = np.mean(signal[frm:to])
    return signal_new

def upsample_signal(signal , up_sample_factor):
    d = 1.0/up_sample_factor
    return np.interp(np.arange(0, len(signal), d), np.arange(0, len(signal)), signal)

def shift_calc(time_signal , t0):
    time_min_diff = np.min(np.abs(time_signal - t0))
    shift = np.where(np.abs(time_signal - t0) == time_min_diff)[0][0]
    return shift
#################################################################################
KVH = 0
TPI = 1
def load_IMU_dataset(imu_file_path, IMU_CHOICE, Freq_INS):
    imu_data = sio.loadmat(imu_file_path)
    #get all important data
    #Time
    time_imu = imu_data.get('IMU_second')[:,0] 
    #get imu frequency
    Freq_imu = np.round(1 / np.mean(time_imu[1:] - time_imu[:-1]))
    print('imu frequency = ' , Freq_imu , ' Hz')
    down_sampleBin_imu = int(np.round(Freq_imu / Freq_INS))
    print('downsampling window  = ' , down_sampleBin_imu)
    #downsample time
    time_imu = down_sample(time_imu , down_sampleBin_imu)
    #-----------------------------------------------
    #accelorometer data
    fx_imu = imu_data['f']['x'][0,0][:,0]
    fy_imu = imu_data['f']['y'][0,0][:,0] 
    if IMU_CHOICE == KVH:
        fy_imu = fy_imu* -1 # special case for the accelerometer of imu is inverted
    fz_imu = imu_data['f']['z'][0,0][:,0]
    #downsample acc
    fx_imu = down_sample(fx_imu , down_sampleBin_imu)
    fy_imu = down_sample(fy_imu , down_sampleBin_imu)
    fz_imu = down_sample(fz_imu , down_sampleBin_imu)
    #-----------------------------------------------
    #Gyro data
    wx_imu = imu_data['w']['x'][0,0][:,0]
    wy_imu = imu_data['w']['y'][0,0][:,0] 
    if IMU_CHOICE == KVH:
        wy_imu = wy_imu* -1 # special case for the accelerometer of imu is inverted
    wz_imu = imu_data['w']['z'][0,0][:,0]
    
    if IMU_CHOICE == TPI:
        wx_imu = wx_imu - np.mean(wx_imu[200:800])
        wy_imu = wy_imu - np.mean(wy_imu[200:800])
        wz_imu = wz_imu - np.mean(wz_imu[200:800])
        
        fx_imu = fx_imu - np.mean(wx_imu[200:800])
        fy_imu = fy_imu - np.mean(wy_imu[200:800])
        fz_imu = fz_imu - np.mean(wz_imu[200:800])
    
    
    #downsample acc
    wx_imu = down_sample(wx_imu , down_sampleBin_imu)
    wy_imu = down_sample(wy_imu , down_sampleBin_imu)
    wz_imu = down_sample(wz_imu , down_sampleBin_imu)
    
    imu_data = dict()
    imu_data['time'] = time_imu
    imu_data['fx'] = fx_imu
    imu_data['fy'] = fy_imu
    imu_data['fz'] = fz_imu
    imu_data['wx'] = wx_imu
    imu_data['wy'] = wy_imu
    imu_data['wz'] = wz_imu   
    
    return imu_data
#################################################################################
######################### Loading Ref Data ######################################
#################################################################################
def load_REF_dataset(ref_file_path, Freq_INS):
    ref_data = sio.loadmat(ref_file_path)
    #get all important data
    #Time
    time_ref = ref_data.get('INS_second')[:,0] 
    #----------------------------------------------------------
    # position info
    long_ref = ref_data.get('INS_Long')[:,0] * np.pi / 180.0
    lat_ref = ref_data.get('INS_Lat')[:,0] * np.pi / 180.0
    alt_ref = ref_data.get('INS_Alt')[:,0] 
    #----------------------------------------------------------
    #Attitude info
    roll_ref = ref_data.get('INS_Roll')[:,0] * np.pi / 180.0
    pitch_ref = ref_data.get('INS_Pitch')[:,0] * np.pi / 180.0
    azimuth_ref = ref_data.get('INS_Azi')[:,0] * np.pi / 180.0
    #----------------------------------------------------------
    #Speed Info
    Ve_ref = ref_data.get('INS_ve')[:,0] 
    Vn_ref = ref_data.get('INS_vn')[:,0] 
    Vu_ref = ref_data.get('INS_vu')[:,0] 
    #get ref frequency
    Freq_ref = np.round(1 / np.mean(time_ref[1:] - time_ref[:-1]))
    if Freq_ref > Freq_INS:
        sampling_factor = int(np.round(Freq_ref / Freq_INS))
        sample_fn = down_sample
    elif Freq_ref < Freq_INS:
        sampling_factor = int(Freq_INS / Freq_ref)
        sample_fn = upsample_signal
    else:
        sampling_factor = 1
    print('Ref frequency = ' , Freq_ref , ' Hz')
    print('Sampling Factor = ' , sampling_factor)
    if sampling_factor != 1:
        time_ref = sample_fn(time_ref,sampling_factor)
        #position 
        long_ref = sample_fn(long_ref,sampling_factor)
        lat_ref = sample_fn(lat_ref,sampling_factor)
        alt_ref = sample_fn(alt_ref,sampling_factor)
        #Attitude 
        roll_ref = sample_fn(roll_ref,sampling_factor)
        pitch_ref = sample_fn(pitch_ref,sampling_factor)
        azimuth_ref = sample_fn(azimuth_ref,sampling_factor)
        #Speed Info
        Ve_ref = sample_fn(Ve_ref,sampling_factor)
        Vn_ref = sample_fn(Vn_ref,sampling_factor)
        Vu_ref = sample_fn(Vu_ref,sampling_factor)

    ref_data = {'time':time_ref , 'lat':lat_ref,'lon':long_ref , 'alt':alt_ref , 
                'roll':roll_ref,'pitch':pitch_ref , 'azimuth':azimuth_ref,
                'Ve':Ve_ref, 'Vn':Vn_ref, 'Vu':Vu_ref}
    return ref_data
#################################################################################
######################### Loading Odo Data ######################################
#################################################################################
def load_ODO_dataset(odo_dir_path, speed_file_name ='odo_second.mat' , time_file_name ='CarChip_Speed.mat' ):
    #Odometry data folder
    #Load the time [ 1 HZ]
    time_odo = sio.loadmat(os.path.join(odo_dir_path , speed_file_name  ))
    #convert to a vector
    time_odo = time_odo['odo_second'][0,:]
    #get odo frequency
    Freq_odo = 1 / np.mean(time_odo[1:] - time_odo[:-1])
    print('Odometry frequency = ' , Freq_odo , ' Hz')
    #Load the speed in m/s
    speed_odo = sio.loadmat(os.path.join(odo_dir_path , time_file_name))
    #convert to a vector
    speed_odo = speed_odo['CarChip_Speed'][:,0]
    #calculate acceleration of the odometer
    acc_odo = np.zeros(speed_odo.shape)
    acc_odo[1:] = (speed_odo[1:] - speed_odo[:-1]) * Freq_odo # (v[t+1] - v[t]) / delta_t    
    
    odo_data = {'time':time_odo,'v_odo':speed_odo , 'a_odo':acc_odo}
    return odo_data
#################################################################################
def load_GPS_dataset(gps_file_path):
    #Note this data set first samples needs to be skipped
    skip = 1
    gp_data = sio.loadmat(gps_file_path)
    #get all important data
    gps_data = dict()
    #Time
    time = gp_data.get('GP_second')[skip:,0]     
    #get gp frequency
    Freq = np.round(1 / np.mean(time[1:] - time[:-1]))
    print('GP frequency = ' , Freq , ' Hz')
    gps_data['time'] = time
    #-----------------------------------------------    
    #Position Data
    gps_data['lat'] = gp_data.get('GP_Lat')[skip:,0] * np.pi / 180.0
    gps_data['lat_std'] = gp_data.get('GP_Lat_std')[skip:,0] 
    
    gps_data['lon'] = gp_data.get('GP_Long')[skip:,0] * np.pi / 180.0
    gps_data['lon_std'] = gp_data.get('GP_Long_std')[skip:,0] 
    
    gps_data['alt'] = gp_data.get('GP_Alt')[skip:,0] 
    gps_data['alt_std'] = gp_data.get('GP_Alt_std')[skip:,0] 
    #-----------------------------------------------    
    #Velocity Data
    gps_data['Ve'] = gp_data.get('GV_ve')[skip:,0]     
    gps_data['Vn'] = gp_data.get('GV_vn')[skip:,0] 
    gps_data['Vu'] = gp_data.get('GV_vu')[skip:,0] 
    
    return gps_data
#################################################################################
def load_BGPS_dataset(gps_pos_file_path , gps_vel_file_path):
    gp_data = sio.loadmat(gps_pos_file_path)
    #get all important data
    gps_data = dict()
    #Time
    time = gp_data.get('BGP_second')[:,0]     
    #get gp frequency
    Freq = np.round(1 / np.mean(time[1:] - time[:-1]))
    print('BGP frequency = ' , Freq , ' Hz')
    gps_data['time'] = time
    #-----------------------------------------------    
    #Position Data
    gps_data['lat'] = gp_data.get('BGP_Lat')[:,0] * np.pi / 180.0
    gps_data['lat_std'] = gp_data.get('BGP_Lat_std')[:,0] 
    
    gps_data['lon'] = gp_data.get('BGP_Long')[:,0] * np.pi / 180.0
    gps_data['lon_std'] = gp_data.get('BGP_Long_std')[:,0] 
    
    gps_data['alt'] = gp_data.get('BGP_Alt')[:,0] 
    gps_data['alt_std'] = gp_data.get('BGP_Alt_std')[:,0] 
    #-----------------------------------------------   
    gp_data = sio.loadmat(gps_vel_file_path)
    #Velocity Data
    gps_data['Ve'] = gp_data.get('BGV_ve')[:,0]
    gps_data['Vn'] = gp_data.get('BGV_vn')[:,0]
    gps_data['Vu'] = gp_data.get('BGV_vu')[:,0]
    return gps_data
#################################################################################
def load_Raw_GPS_dataset(gps_file_path):
    #Note this data set first samples needs to be skipped
    skip = 1
    gp_data = sio.loadmat(gps_file_path)
    #get all important data
    gps_data = dict()
    #Time
    time = gp_data.get('GP_second')[skip:,0]     
    #get gp frequency
    Freq = np.round(1 / np.mean(time[1:] - time[:-1]))
    print('GP frequency = ' , Freq , ' Hz')
    gps_data['time'] = time
    #-----------------------------------------------    
    #Position Data
    gps_data['lat'] = gp_data.get('GP_Lat')[skip:,0] * np.pi / 180.0
    gps_data['lat_std'] = gp_data.get('GP_Lat_std')[skip:,0] 
    
    gps_data['lon'] = gp_data.get('GP_Long')[skip:,0] * np.pi / 180.0
    gps_data['lon_std'] = gp_data.get('GP_Long_std')[skip:,0] 
    
    gps_data['alt'] = gp_data.get('GP_Alt')[skip:,0] 
    gps_data['alt_std'] = gp_data.get('GP_Alt_std')[skip:,0] 
    #-----------------------------------------------    
    #Velocity Data
    gps_data['Ve'] = gp_data.get('GV_ve')[skip:,0]     
    gps_data['Vn'] = gp_data.get('GV_vn')[skip:,0] 
    gps_data['Vu'] = gp_data.get('GV_vu')[skip:,0] 
    
    return gps_data
################################################################################
#Data Skipping
def crop_data(data, skip_time,run_time):
    dt = data['time'][1] - data['time'][0]
    F = int(np.round(1/dt))
    if skip_time != None:
        #find index of this time
        i_start = shift_calc(data['time'] , skip_time)
        print('skipping till = ',skip_time , ' = ' , i_start , ' samples')
    else:        
        i_start = 0
    #-----------------------------------   
    #Skipping here
    if run_time != None:        
        samples_to = i_start + int(run_time *60*F)
        print('Running for till = ',samples_to , 'samples')
    else:
        samples_to = len(data['time']) 
        
        
    #apply shift
    for key in data:
        data[key] = data[key][i_start:samples_to]    

    return data


#################################################################################
######################### Data Syncronization ###################################
#################################################################################
def synchronize_general(data,time_from,time_to):
    dt = data['time'][1] - data['time'][0]
    if (time_from - data['time'][0]) > dt:
        #calculate shift
        _start_shift = shift_calc(data['time'] , time_from)
                    
        #apply shift
        for key in data:
            data[key] = data[key][_start_shift:]
        
    if (data['time'][-1] - time_to) > dt :
        #calculate shift
        _end_shift = shift_calc(data['time'] , time_to) + 1
        #apply shift
        for key in data:
            data[key] = data[key][:_end_shift]
    
def syncronize_INS(ref_data , imu_data):
    #--------------------------------------------------
    #perform needed time shifting to synchronize signals
    #The start time should be the latest time in the data time frames
    time_from = max(ref_data['time'][0], imu_data['time'][0])
    #The end time should be the earliest time in the data time frames
    time_to = min(ref_data['time'][-1], imu_data['time'][-1])    
    #shifting all needed signals
    #-----------------------------------------------
    #1. Reference data
    synchronize_general(ref_data,time_from,time_to)
    print('* Ref -> ',len(ref_data['time']) , ref_data['time'][0] , ref_data['time'][-1])
    #------------------------------------------------
    #2. IMU data
    synchronize_general(imu_data,time_from,time_to)
    print('* IMU -> ',len(imu_data['time']) , imu_data['time'][0] , imu_data['time'][-1])

    return ref_data , imu_data
##################################################################################################################################################################
def syncronize_RISS(ref_data , imu_data, odo_data):
    #--------------------------------------------------
    #perform needed time shifting to synchronize signals
    #The start time should be the latest time in the data time frames
    time_from = max(ref_data['time'][0], imu_data['time'][0], odo_data['time'][0])
    #The end time should be the earliest time in the data time frames
    time_to = min(ref_data['time'][-1], imu_data['time'][-1], odo_data['time'][-1])
    #shifting all needed signals
    #-----------------------------------------------
    #1. Reference data
    synchronize_general(ref_data,time_from,time_to)
    print('* REF -> ',len(ref_data['time']) , ref_data['time'][0] , ref_data['time'][-1])
    #------------------------------------------------
    #2. IMU data
    synchronize_general(imu_data,time_from,time_to)
    print('* IMU -> ',len(imu_data['time']) , imu_data['time'][0] , imu_data['time'][-1])
    #------------------------------------------------
    #2. Odo data
    synchronize_general(odo_data,time_from,time_to)
    print('* ODO -> ',len(odo_data['time']) , odo_data['time'][0] , odo_data['time'][-1])

    return ref_data , imu_data,odo_data
##################################################################################################################################################################
def syncronize_RISS_GPS(ref_data , imu_data, odo_data, gps_data):
    #--------------------------------------------------
    #perform needed time shifting to synchronize signals
    #The start time should be the latest time in the data time frames
    time_from = max(ref_data['time'][0], imu_data['time'][0], odo_data['time'][0],gps_data['time'][0])
    #The end time should be the earliest time in the data time frames
    time_to = min(ref_data['time'][-1], imu_data['time'][-1], odo_data['time'][-1],gps_data['time'][-1])
    #shifting all needed signals
    #-----------------------------------------------
    #1. Reference data
    synchronize_general(ref_data,time_from,time_to)
    print('* REF -> ',len(ref_data['time']) , ref_data['time'][0] , ref_data['time'][-1])
    #------------------------------------------------
    #2. IMU data
    synchronize_general(imu_data,time_from,time_to)
    print('* IMU -> ',len(imu_data['time']) , imu_data['time'][0] , imu_data['time'][-1])
    #------------------------------------------------
    #3. Odo data
    synchronize_general(odo_data,time_from,time_to)
    print('* ODO -> ',len(odo_data['time']) , odo_data['time'][0] , odo_data['time'][-1])
    #------------------------------------------------
    #4. GPS data
    synchronize_general(gps_data,time_from,time_to)
    print('* GPS -> ',len(gps_data['time']) , gps_data['time'][0] , gps_data['time'][-1])

    return ref_data , imu_data,odo_data,gps_data
##################################################################################################################################################################
def syncronize_INS_GPS(ref_data , imu_data,gps_data):
    #--------------------------------------------------
    #perform needed time shifting to synchronize signals
    #The start time should be the latest time in the data time frames
    time_from = max(ref_data['time'][0], imu_data['time'][0],gps_data['time'][0])
    #The end time should be the earliest time in the data time frames
    time_to = min(ref_data['time'][-1], imu_data['time'][-1],gps_data['time'][-1])    
    #shifting all needed signals
    #-----------------------------------------------
    #1. Reference data
    synchronize_general(ref_data,time_from,time_to)
    print('* Ref -> ',len(ref_data['time']) , ref_data['time'][0] , ref_data['time'][-1])
    #------------------------------------------------
    #2. IMU data
    synchronize_general(imu_data,time_from,time_to)
    print('* IMU -> ',len(imu_data['time']) , imu_data['time'][0] , imu_data['time'][-1])
    #------------------------------------------------
    #3. GPS data
    synchronize_general(gps_data,time_from,time_to)
    print('* GPS -> ',len(gps_data['time']) , gps_data['time'][0] , gps_data['time'][-1])

    return ref_data , imu_data , gps_data

##################################################################################################
###################################### Testing Section ###########################################
##################################################################################################
import pandas as pd
def prepare_data_for_ML(Folder_Path,IMU_CHOICE=KVH):
    #Reference Data loading
    ref_file_path = r'E:\Eng Work\6-PHD Study\Second Semster\EEE523\data set\NovAtel\3.INSPVA_Reference.mat'
    ref_data = load_REF_dataset(ref_file_path, Freq_INS=5)
    #-------------------------------------------------------
    # IMU Data Loading
    #KVH
    imu_file_path = r'E:\Eng Work\6-PHD Study\Second Semster\EEE523\data set\NovAtel\1.RAWIMU.mat'
    imu_data = load_IMU_dataset(imu_file_path, IMU_CHOICE=KVH, Freq_INS=20)
    #TPI
    #imu_file_path = r'E:\Eng Work\6-PHD Study\Second Semster\EEE523\data set\TPI\TPI_data_interpolated2.mat'
    #imu_data = load_IMU_dataset(imu_file_path, IMU_CHOICE=TPI, Freq_INS=1)
    #-------------------------------------------------
    # ODO Data Loading
    odo_dir_path = r'E:\Eng Work\6-PHD Study\Second Semster\EEE523\data set\OBDII_data'
    odo_data = load_ODO_dataset(odo_dir_path)
    #---------------------------------------------------------
    # GPS Data Loading
    #Novatel
    gps_dir_path = r"E:\Eng Work\6-PHD Study\Second Semster\EEE523\data set\NovAtel"
    gps_data = load_BGPS_dataset(os.path.join(gps_dir_path,"2.BESTGPSPOS.mat"),os.path.join(gps_dir_path,"2.BESTGPSVEL.mat"))
    
    
    ref_data , imu_data, odo_data, gps_data = syncronize_RISS_GPS(ref_data , imu_data, odo_data, gps_data)
    
    #Save reference
    df = pd.DataFrame.from_dict(ref_data)
    df.to_csv(os.path.join(Folder_Path,'reference_5Hz.csv'),index=False)
    df = pd.DataFrame.from_dict(imu_data)
    df.to_csv(os.path.join(Folder_Path,'imu_KVH_20Hz.csv'),index=False)
    df = pd.DataFrame.from_dict(odo_data)
    df.to_csv(os.path.join(Folder_Path,'odo.csv'),index=False)
    df = pd.DataFrame.from_dict(gps_data)
    df.to_csv(os.path.join(Folder_Path,'gps_Novatel.csv'),index=False)
################################################################################
if __name__ == "__main__":
    # Testing Section    
    #Folder_Path = r'E:\Eng Work\6-PHD Study\Second Semster\ML Course\Project\DataSet\KVH'
    #Folder_Path = r"E:\OneDrive - Queen's University\ELEC825\DataSet Kingston\KVH"
    #prepare_data_for_ML(Folder_Path)
    
    #Testing GPS Raw Data Loading
    Folder_Path = r"E:\Eng Work\6-PHD Study\Second Semster\EEE523\Project_GPS"
    file_path = os.path.join(Folder_Path,"GPS_data_combined_Nova.mat")
    
    