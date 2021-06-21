from IntializationParameters import InitINS
import numpy as np
import math

class Mechanization (InitINS): 
    def __init__ (self, latitude, longitude, altitude, roll, pitch, yaw):
        super().__init__(latitude, longitude, altitude, roll, pitch, yaw)
        self.Rb_l = self.Init_Rbl
        self.Quatrenion = self.Init_Quatrenion
        self.Localg = self.Init_Localg

    #### Calculate Radii of curvature
    def RadiiM (self, latitude):
        self.Rm =  ( self.a *(1-self.e2) ) / ( (1- self.e2 * np.sin(latitude) * np.sin(latitude) )**(1.5) )
    def RadiiN (self,latitude):
        self.Rn = self.a / np.sqrt ( 1- ( self.e2 * np.sin (latitude) * np.sin (latitude) ) )
    ### Calculate Transformation from ECF to Local level frame
    def R_EL (self, latitude, longitude) :
        """
        Input:  Latitude, and Longitude
        Output: Updates the rotation matrix from the ECEF to Local level frame
        """
        # Rotation from ECEF to local level frame
        ### Note: 1- The numpy trigonometric functions take radians as an input.
        ##### First Row  #####
        self.rel11 = ( - np.sin (0) * np.sin(latitude) * np.cos (longitude) ) -  ( np.cos (0) * np.sin (longitude) ) 

        self.rel12 = (-np.sin(0) * np.sin (latitude) * np.sin (longitude) ) + ( np.cos(0) * np.cos (longitude) )

        self.rel13 = np.sin(0) * np.cos (latitude)
        ##### Second Row #####
        self.rel21 = ( - np.cos(0) * np.sin(latitude) * np.cos (longitude) ) + ( np.sin(0) * np.sin(longitude) )

        self.rel22 = ( -np.cos(0) * np.sin (latitude) * np.sin (longitude) ) - ( np.sin(0) * np.cos(longitude) )

        self.rel23 = np.cos (0) * np.cos(latitude)
        ##### Third Row ######
        self.rel31 = ( np.cos (latitude) * np.cos (longitude) )

        self.rel32 =  np.cos (latitude) * np.sin (longitude)

        self.rel33 = np.sin (latitude)
        ##### Matrix formulation ######
        self.Re_l = np.array ([[self.rel11,self.rel12,self.rel13],[self.rel21,self.rel22,self.rel23],[self.rel31,self.rel32,self.rel33]])

        
    def WIE_L (self) :
        '''
        Input: None
        Output: transforms the earth angular velocity into the local level frame using the R_EL
        '''
        we_transform = np.transpose(np.array([0,0,self.We]))
        self.wie_l = np.dot(self.Re_l, we_transform)

    def WEL_L (self, latitude, altitude):
        '''
        Input: Laitiude, Altitude
        Output: calculates the angular velocities vector using the velocity vector
        '''
        wel_l1 = - self.vl [1] / (self.Rm + altitude)
        wel_l2 = self.vl [0]/ ( self.Rn + altitude)
        wel_l3 = ( self.vl [0] * np.tan (latitude) )  /( self.Rn + altitude)
        self.wel_l = np.transpose ( np.array([wel_l1, wel_l2, wel_l3]))
    def WLB_B (self, wx, wy, wz) :
        wib_b = np.transpose(np.array([wx,wy,wz]))
        self.wlb_b = wib_b -  np.dot( np.transpose(self.Rb_l) , ( self.wel_l + self.wie_l ) )
    
    def SkewMatrix_WLB_B (self):
        '''
        Input: None
        Output: The skewed Wlb_b matrix
        It calculates the skew matrix of the angular velocity components
        '''
        self.skewmatrix_wlb_b = np.array([[0,self.wlb_b [2] , -self.wlb_b[1] , self.wlb_b[0]]
                                        ,[-self.wlb_b[2], 0, self.wlb_b[0], self.wlb_b[1]]
                                        ,[self.wlb_b[2], - self.wlb_b[0], 0 ,self.wlb_b[2]]
                                        ,[-self.wlb_b[0], - self.wlb_b[1], -self.wlb_b[2], 0]])
 
    def UpdateQuatrenion (self, delta_time):
        '''
        Input: The change in time (Delta Time)
        Output: Updates the Quatrenion vector
        '''
        Q_dot = 0.5 * np.dot (self.skewmatrix_wlb_b, self.Quatrenion)
        self.Quatrenion += (delta_time * Q_dot)
        self.Quatrenion /= np.norm(self.Quatrenion)

    def UpdateRBL (self):
        '''
        Input: None
        Ouput: Updates the Rotation from body to local level frame using the updated Quatrnions
        '''
         ##### First Row  #####
        self.rbl11 = ( self.Quatrenion[0]**2 ) - ( self.Quatrenion[1]**2 ) - ( self.Quatrenion[2] * ( self.Quatrenion[3]**2 ) )

        self.rbl12 = 2 * ( ( self.Quatrenion[0] * self.Quatrenion[1] ) - ( self.Quatrenion[2] *  self.Quatrenion[3] ) )

        self.rbl13 = 2 * ( ( self.Quatrenion[0] * self.Quatrenion[2] ) + ( self.Quatrenion[1] * self.Quatrenion[3] ) ) 

        ##### Second Row #####
        self.rbl21 =  2 * ( ( self.Quatrenion[0] * self.Quatrenion[1] ) + ( self.Quatrenion[2] *  self.Quatrenion[3] ) )

        self.rbl22 =  - (self.Quatrenion[0]**2) + (self.Quatrenion[1]**2 ) - (self.Quatrenion[2]**2) + ( self.Quatrenion[3]**2 )

        self.rbl23 = 2 * ( ( self.Quatrenion[1] * self.Quatrenion[2] ) - ( self.Quatrenion[0] * self.Quatrenion[3] ) )

        ##### Third Row ######
        self.rbl31 = 2 * ( ( self.Quatrenion[0] * self.Quatrenion[2] ) - ( self.Quatrenion[1] * self.Quatrenion[3] ) )

        self.rbl32 = 2 * ( ( self.Quatrenion[1] * self.Quatrenion[3] ) + ( self.Quatrenion[0] * self.Quatrenion[4] ) )

        self.rbl33 = - (self.Quatrenion[0]**2) - (self.Quatrenion[1]**2 ) + (self.Quatrenion[2]**2) + ( self.Quatrenion[3]**2 )

        ##### Matrix formulation ######
        self.Rb_l = np.array ([[self.rbl11,self.rbl12,self.rbl13],[self.rbl21,self.rbl22,self.rbl23],[self.rbl31,self.rbl32,self.rbl33]])

    def UpdateAttitude (self):
        '''
        Input: None
        Output: Returns the updated attitude components (roll, pitch yaw (azimuth)) 
        '''
        self.pitch  = np.asin (self.rbl32)
        self.roll   = np.atan2(- self.rbl31 , self.rbl33 )
        self.yaw    = np.atan2(- self.rbl12 , self.rbl22 )

        if self.yaw < 0:
            self.yaw += (2*np.pi)
        elif self.yaw >= (2*np.pi):
            self.yaw += (- 2*np.pi)
        
        self.azimuth = (2*np.pi) - self.yaw

    def OMEGA_IE_L (self):
        self.omega_ie_l = np.array([[0, - self.wie_l[2], self.wie_l[1]]
                                ,[self.wie_l[2], 0 , - self.wie_l[0]]
                                ,[- self.wie_l[1], self.wie_l[0], 0]])
    
    def OMEGA_EL_L (self):
        self.omega_el_l = np.array([[0, - self.wel_l[2], self.wel_l[1]]
                                ,[self.wel_l[2], 0 , - self.wel_l[0]]
                                ,[- self.wel_l[1], self.wel_l[0], 0]])
    def UpdateAccelerometers (self, fx, fy, fz):
        self.fb = np.transpose(np.array([fx,fy,fz]))

    def UpdateG (self,g):
        self.gVector = np.transpose(np.array([ 0, 0, -g])) 
    def UpdateDeltaVelocity (self, delta_time):
        self.Delta_Vl_current = ( np.dot(self.Rb_l,self.fb) *  delta_time ) - ( np.dot( (2 * self.omega_ie_l + self.omega_el_l), self.vl) * delta_time ) + (self.gVector * delta_time)
    
    def UpdateVelocity (self):
        self.vl =+ (0.5 * (self.Delta_Vl_current + self.Delta_Vl_previous ))
        self.Delta_Vl_previous = self.Delta_Vl_current

    def UpdatePosition (self, prev_latitude, prev_longitude, prev_altitude, ve_prev, vn_prev, vu_prev, delta_time):
        ve, vn, vu = self.vl[0], self.vl[1], self.vl[2]
        self.longitude = prev_longitude + ( 0.5 * ( ( (ve + ve_prev) * delta_time ) / ((self.Rn + prev_altitude) * np.cos(prev_latitude))))
        self.latitude = prev_latitude + 0.5 * ( ( vn + vn_prev) * delta_time ) / (self.Rm + prev_altitude) 
        self.altitude = prev_altitude + 0.5 * ((vu + vu_prev) * delta_time)


    def compile (self, dt, wx, wy, wz
                , fx, fy,fz 
                ,latitude,longitude,altitude,
                 ve_Prev, vn_Prev, vu_Prev,
                 g ):
        '''
        Input: Sample Time (dt),
        Current Gyroscopes components (wx,wy,wz),
        Current accelerometer components (fx, fy ,fz),
        Previous epoch's position components (latitude, longitude, altitude),
        Previous epoch's velocity components (ve_prev, vn_prev, vu_prev),
        The gravity acceleration (g).

        Output: Returns the current epoch's (position components, velocity components, attitude components)

        This method performs the complete INS mechanization in a single epoch
        '''
        ###### Stage 1: updating the earth's radii of Curvature  #####
        self.RadiiM(latitude)
        self.RadiiN(latitude)

        ##### Stage 2: updating the angular velocities #########
        self.R_EL (latitude, longitude)
        self.WIE_L()
        self.WEL_L(latitude,altitude)
        self.WLB_B(wx,wy,wz)
        self.SkewMatrix_WLB_B()
        #### Stage 3: updating quatrenions ##################
        self.UpdateQuatrenion(dt)
        #### Stage 4: updating attitude components ########
        self.UpdateRBL()
        self.UpdateAttitude()
        ##### Stage 5: updating velocity components #######
        self.OMEGA_IE_L()
        self.OMEGA_EL_L()
        self.UpdateAccelerometers(fx,fy,fz)
        self.UpdateG(g)
        self.UpdateDeltaVelocity (dt)
        self.UpdateVelocity ()
        ##### Stage 6: updating position components #######
        self.UpdatePosition (latitude, longitude, altitude, ve_Prev, vn_Prev, vu_Prev, dt)

        ##### Return the updated full navigation components #####
        return self.latitude, self.longitude, self.altitude , self.vl[0],self.vl[1], self.vl[2], self.roll, self.pitch, self.azimuth,self.yaw
                



       
