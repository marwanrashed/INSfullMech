import numpy as np 
from numpy.core.fromnumeric import shape

'''
Property of: NavINST Laboratory
Author: Marwan A. Rashed
Description: The intialization parameters for the full INS mechanization
'''
class InitINS (): 
    def __init__ (self, latitude, longitude, altitude, roll, pitch, azimuth):
        '''
        Input : the intital 3D position components (Latitude, Longitude, Altitude)
                3D Attitude Angles (Roll, Pitch, Yaw)
        Output: None, Initiate local intital variables to intiate the mechanization 
        '''
        # Variables declaration
        ## Position components
        self.Init_Lat, self.Init_Long, self.Init_Alt =  np.deg2rad(latitude) , np.deg2rad(longitude), altitude
        ## Attitude angles
        self.Init_Roll, self.Init_Pitch, self.Init_Azimuth =   np.deg2rad (roll) , np.deg2rad (pitch) , np.deg2rad(azimuth)
        

        ## INS Mechanization constant variables
        self.a= 6378137 # the constant for Earth Radii Rm and Rn
        self.fa= 1/298.257223563 
        self.b= self.a * (1-self.fa) 
        self.e2 = 1-(self.b**2) / (self.a**2)         
        self.We= 7.292115e-5 # Earth rotation rate rad/sec

        ## Local gravity computation constants
        self.a1, self.a2, self.a3 = 9.7803267714 , 0.0052790414 , 0.0000232718
        self.a4, self.a5, self.a6 = -0.000003087691089 , 0.000000004397731 , 0.000000000000721

        #### Intialization of change in local velocity vector #####
        self._Delta_Vl_previous = np.array ( np.zeros ((3))).transpose() 

    def Init_Rbl (self): 
        """
        Input: None, takes attitude components from the class initialization
        Output: Returns the initial rotation matrix from the body frame to the local level frame
        """
        # Rotation from body frame to local level frame
        ### Note: 1- The numpy trigonometric functions take radians as an input.
        ### 2- These matrix is intiated based on pre-saved entries / not real time. If real time the structure will have to change
        ##### First Row  #####
        self.rbl11 = ( np.cos(self.Init_Roll) * np.cos(self.Init_Azimuth) ) +  ( np.sin(self.Init_Roll) * np.sin(self.Init_Pitch) * np.sin(self.Init_Azimuth) ) 

        self.rbl12 = np.sin(self.Init_Azimuth) * np.cos(self.Init_Pitch)

        self.rbl13 = np.cos(self.Init_Azimuth)*np.sin(self.Init_Roll) - np.sin(self.Init_Azimuth)*np.sin(self.Init_Pitch)*np.cos(self.Init_Roll)
        ##### Second Row #####
        self.rbl21 = - np.sin(self.Init_Azimuth) * np.cos(self.Init_Roll) + np.cos(self.Init_Azimuth)*np.sin(self.Init_Pitch)*np.sin(self.Init_Roll)

        self.rbl22 = np.cos(self.Init_Azimuth)*np.cos(self.Init_Pitch) 

        self.rbl23 = - np.sin(self.Init_Azimuth)*np.sin(self.Init_Roll) - np.cos(self.Init_Azimuth)*np.sin(self.Init_Pitch)*np.cos(self.Init_Roll)
        ##### Third Row ######
        self.rbl31 = - np.cos(self.Init_Pitch)*np.sin(self.Init_Roll) 

        self.rbl32 = np.sin(self.Init_Pitch)

        self.rbl33 = np.cos(self.Init_Roll) * np.cos(self.Init_Pitch)
        ##### Matrix formulation ######
        self._Rb_l = np.ndarray (shape= (3,3), dtype=float ,buffer= np.array([[self.rbl11,self.rbl12,self.rbl13],[self.rbl21,self.rbl22,self.rbl23],[self.rbl31,self.rbl32,self.rbl33]]))
        print ("Intial R_Bl:",self._Rb_l)

    def Init_Quatrenion (self):
        """
        Input: None, takes Rotation Body to local matrix components from the class shared variables
        Output: Returns the intial quatrenion vector
        """
        # Quatrenions Formulation
        ##### Quatrnions declaration #####
        Fourth_Quatrenion = np.sqrt (1 + self.rbl11 + self.rbl22 + self.rbl33) * 0.5
        First_Quatrenion = 0.25 * (self.rbl32 - self.rbl23 ) / Fourth_Quatrenion
        Second_Quatrenion = 0.25*( self.rbl13 - self.rbl31 ) / Fourth_Quatrenion
        Third_Quatrenion = 0.25*( self.rbl21 - self.rbl12) / Fourth_Quatrenion
        #### Quatrenions Vector formulation #####
        self._Quatrenion =  np.array([First_Quatrenion , Second_Quatrenion , Third_Quatrenion , Fourth_Quatrenion]).transpose()
        print ("Intial Quaterenion before correction:", self._Quatrenion)
        self._Quatrenion = self._Quatrenion / np.linalg.norm(self._Quatrenion)
        print ("Intial Quaterenion after correction:",self._Quatrenion)
        # return self._Quatrenion

    def Init_Localg (self):
        '''
        Input: None
        Output: Calculates the local gravity component, which is shared across the object.
        '''
        Proj_Lat = np.sin( self.Init_Lat ) # Projection of the latitude after degree to raddian transformation
        self._Local_g = self.a1  * (1 + (self.a2* (Proj_Lat**2) ) + (   self.a3  * (Proj_Lat**4) ) ) + ( ( self.a4 + (self.a5 * (Proj_Lat**2)) ) *  self.Init_Alt ) + self.a6 * (self.Init_Alt**2)
        # return self._Local_g

    def Init_Velocity (self, Initial_ve = 0 , Initial_vn= 0, Initial_vu = 0):
        '''
        Input: The Initial 3D velocity vector 
        Output: None, updates the local shared velocity vector variable.
        '''
        self._vl = np.array([Initial_ve, Initial_vn, Initial_vu]).transpose()
        print ("Intial Velocities", self._vl )




        









