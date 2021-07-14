#################################################################################
############################ Plotting Funs ######################################
#################################################################################
import matplotlib
import numpy as np
from matplotlib import pyplot as plt
matplotlib.rcParams.update({'font.size': 16})

def plot_to_compare(sig_ref,time_ref, sig_imu,time_imu, title_text , ylabel , solution_label , gps_outage_info = None ):
        fig = plt.figure(dpi=100 , figsize=(7,4))
        time_ref_plotting = (time_ref - time_ref.min())/60.0
        time_imu_plotting = (time_imu - time_imu.min())/60.0
        
        plt.plot(time_imu_plotting,sig_imu,'b',label = solution_label)
        plt.plot(time_ref_plotting,sig_ref,'r--',label = 'Ref')
        plt.title(title_text)
        plt.xlabel('time [mins]')
        plt.ylabel(ylabel)
        plt.legend()
        plt.grid()
        #plt.ticklabel_format(useOffset=False)
        
        if gps_outage_info != None:
            for s , e in zip(gps_outage_info['outage_stime'],gps_outage_info['outage_etime']):
                plt.axvspan(s, e, color='yellow', alpha=0.8)
        
        plt.tight_layout()
        plt.show()
    
def plot_error(sig_ref,time_ref, sig_imu,time_imu , title_text , ylabel, gps_outage_info = None):
    fig = plt.figure(dpi=100 , figsize=(7,4))
    time_imu_plotting = (time_imu - time_imu.min())/60.0
    abs_err = np.abs(sig_ref - sig_imu)
    plt.plot(time_imu_plotting,abs_err,'k--')
    plt.title(title_text)
    plt.xlabel('time [mins]')
    plt.ylabel(ylabel)
    plt.grid()
    if gps_outage_info != None:
            for s , e in zip(gps_outage_info['outage_stime'],gps_outage_info['outage_etime']):
                plt.axvspan(s, e, color='yellow', alpha=0.8)
    plt.tight_layout()
    plt.show()
    
def plot_distance( distance_error,time_imu , title_text, gps_outage_info = None):
    fig = plt.figure(dpi=100 , figsize=(7,4))
    plt.title(title_text)
    time_imu_plotting = (time_imu - time_imu.min())/60.0
    plt.plot(time_imu_plotting, distance_error, 'k--')
    plt.ylabel('[m]')
    plt.xlabel('time [mins]')
    plt.grid()
    if gps_outage_info != None:
            for s , e in zip(gps_outage_info['outage_stime'],gps_outage_info['outage_etime']):
                plt.axvspan(s, e, color='yellow', alpha=0.8)
    plt.tight_layout()
    plt.show()

def compute_error (yp,y):
    return np.abs(yp-y)
    
def compute_RMSE(yp , y):
    return np.sqrt(np.mean((yp - y)**2))

def compute_MaxAE(yp , y):
    return (np.max(np.abs(yp - y)))


def plot_bi_trajectory(sig_ref, sig_imu, title_text, solution_label):
        fig = plt.figure(dpi=100 , figsize=(7,4))
        
        plt.plot(sig_imu['lon'],sig_imu['lat'],'b',label = solution_label)
        plt.plot(sig_ref['lon'],sig_ref['lat'],'r--',label = 'Ref')
        plt.title(title_text)
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        plt.legend()
        plt.grid()
        #plt.ticklabel_format(useOffset=False)

        plt.tight_layout()
        plt.show()

def plot_tri_trajectory(sig_ref, sig_imu, sig_aux, title_text, solution_label):
        fig = plt.figure(dpi=100 , figsize=(7,4))

        plt.plot(sig_aux['lon'],sig_aux['lat'],'b',label = solution_label)
        plt.plot(sig_imu['lon'],sig_imu['lat'],'b',label = solution_label)
        plt.plot(sig_ref['lon'],sig_ref['lat'],'r--',label = 'Ref')
        plt.title(title_text)
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        plt.legend()
        plt.grid()
        #plt.ticklabel_format(useOffset=False)

        plt.tight_layout()
        plt.show()

def full_diagnosis ( ref_data,  sig_data, label):
    solution_label = label
    ##########################Plot Trajectory ########################
    plot_bi_trajectory(ref_data, sig_data, 'Trajectory', solution_label)
    ########################## Plot Position errors ###################
    plot_error(ref_data['lat'],ref_data['time'], sig_data['lat'] ,sig_data['time'], 'Latitude errors', 'error in degrees',  )
    plot_error(ref_data['lon'],ref_data['time'], sig_data['lon'],sig_data['time'], 'Longitude errors', 'error in degrees')
    plot_error(ref_data['alt'],ref_data['time'], sig_data['alt'],sig_data['time'], 'Altitude errors', 'error in meters')
    ########################## Plot Attitude errors ###################
    plot_error(ref_data['roll'],ref_data['time'], sig_data['roll'],sig_data['time'], 'Roll errors', 'error in degrees')
    plot_error(ref_data['pitch'],ref_data['time'], sig_data['pitch'],sig_data['time'], 'Pitch errors', 'error in degrees')
    plot_error(ref_data['azimuth'],ref_data['time'], sig_data['azimuth'],sig_data['time'], 'Azimuth errors', 'error in degree')
    ########################## Plot Velocity errors ###################
    plot_error(ref_data['ve'],ref_data['time'], sig_data['ve'],sig_data['time'], 'East Velocity errors', 'error in m/s')
    plot_error(ref_data['vn'],ref_data['time'], sig_data['vn'],sig_data['time'], 'North Velocity errors', 'error in m/s')
    plot_error(ref_data['vu'],ref_data['time'], sig_data['vu'],sig_data['time'], 'Up Velocity errors', 'error in m/s')
    ######################### Plot in parallel ########################
    ########################## Plot Position  ###################
    plot_to_compare(ref_data['lat'],ref_data['time'], sig_data['lat'],sig_data['time'], 'Latitude ', 'degrees', solution_label)
    plot_to_compare(ref_data['lon'],ref_data['time'], sig_data['lon'],sig_data['time'], 'Longitude ', ' degrees', solution_label)
    plot_to_compare(ref_data['alt'],ref_data['time'], sig_data['alt'],sig_data['time'], 'Altitude', ' meters', solution_label)
    ########################## Plot Attitude errors ###################
    plot_to_compare(ref_data['roll'],ref_data['time'], sig_data['roll'],sig_data['time'], 'Roll ', 'degrees', solution_label)
    plot_to_compare(ref_data['pitch'],ref_data['time'], sig_data['pitch'],sig_data['time'], 'Pitch ', 'degrees', solution_label)
    plot_to_compare(ref_data['azimuth'],ref_data['time'], sig_data['azimuth'],sig_data['time'], 'Azimuth ', ' degree', solution_label)
    ########################## Plot Velocity errors ###################
    plot_to_compare(ref_data['ve'],ref_data['time'], sig_data['ve'],sig_data['time'], 'East Velocity ', 'm/s', solution_label)
    plot_to_compare(ref_data['vn'],ref_data['time'], sig_data['vn'],sig_data['time'], 'North Velocity ', 'm/s', solution_label)
    plot_to_compare(ref_data['vu'],ref_data['time'], sig_data['vu'],sig_data['time'], 'Up Velocity ', 'm/s', solution_label)
  
