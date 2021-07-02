from IntializationParameters import InitINS
import numpy as np
import math

class Mechanization (InitINS): 
    def __init__ (self, latitude, longitude, altitude, roll, pitch, azimuth):
        super().__init__(latitude, longitude, altitude, roll, pitch, azimuth)
        self.Init_Rbl ()
        self.Init_Quatrenion ()
        self.Init_Localg ()

    #### Calculate Radii of curvature
    def RadiiM (self, latitude):
        latitude = np.deg2rad(latitude)
        # print ("Latitude in radian", latitude )
        self._Rm =  ( self.a *(1-self.e2) ) / ( (1- self.e2 * np.sin(latitude) * np.sin(latitude) )**(1.5) )
        # print ("Meredian Radius", self._Rm )
    def RadiiN (self,latitude):
        latitude = np.deg2rad(latitude)
        self._Rn = self.a / np.sqrt ( 1- ( self.e2 * np.sin (latitude) * np.sin (latitude) ) )
        # print ("nORMAL Radius", self._Rn )
    ### Calculate Transformation from ECF to Local level frame
    def R_EL (self, latitude, longitude) :
        """
        Input:  Latitude, and Longitude
        Output: Updates the rotation matrix from the ECEF to Local level frame
        """
        # Rotation from ECEF to local level frame
        ### Note: 1- The numpy trigonometric functions take radians as an input.
        ##### First Row  #####
        latitude  = np.deg2rad(latitude)
        longitude = np.deg2rad (longitude)
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
        self._Re_l = np.ndarray (shape= (3,3), dtype=float
                     ,buffer= np.array ([[self.rel11,self.rel12,self.rel13]
                     ,[self.rel21,self.rel22,self.rel23]
                     ,[self.rel31,self.rel32,self.rel33]]))

        
    def WIE_L (self, latitude) :
        '''
        Input: None
        Output: transforms the earth angular velocity w_e into the local level frame to be the vector WIE_L using the R_EL transformation matrix
        '''
        latitude = np.deg2rad(latitude)
        # we_transform = np.vstack([0,0,self.We])
        # self._wie_l = self._Re_l @ we_transform
        self._wie_l = np.array([0,self.We * np.cos(latitude),self.We * np.sin(latitude)]).transpose()
        # print ("Transformed earth rotation vector",self._wie_l)

    def WEL_L (self, latitude, altitude):
        '''
        Input: Laitiude, Altitude
        Output: calculates the angular velocities vector using the velocity vector
        '''
        latitude = np.deg2rad(latitude)
        wel_l1 = - self._vl [1] / (self._Rm + altitude)
        wel_l2 = self._vl [0]/ ( self._Rn + altitude)
        wel_l3 = ( self._vl [0] * np.tan (latitude) )  /( self._Rn + altitude)
        self._wel_l = np.array ([wel_l1, wel_l2, wel_l3]).transpose()
        # print ("wel_l",self._wel_l)
    def WLB_B (self, wx, wy, wz) :
        wib_b = np.array([wx,wy,wz]).transpose()
        # print ("gyro rates",wib_b)
        self._wlb_b = wib_b -  np.matmul ( np.transpose(self._Rb_l) , ( self._wel_l + self._wie_l ) )
        # print ("WLB_B",self._wlb_b)
    
    def SkewMatrix_WLB_B (self):
        '''
        Input: None
        Output: The skewed Wlb_b matrix
        It calculates the skew matrix of the angular velocity components
        '''
        self._skewmatrix_wlb_b = np.ndarray (shape= (4,4), dtype=float 
                                ,buffer= np.array([[0,self._wlb_b [2] , -self._wlb_b[1] , self._wlb_b[0]]
                                        ,[-self._wlb_b[2], 0, self._wlb_b[0], self._wlb_b[1]]
                                        ,[self._wlb_b[1], - self._wlb_b[0], 0 ,self._wlb_b[2]]
                                        ,[-self._wlb_b[0], - self._wlb_b[1], -self._wlb_b[2], 0]]) )
        # print ("SkewMatrix_WLB_B",self._skewmatrix_wlb_b)
 
    def UpdateQuatrenion (self, delta_time):
        '''
        Input: The change in time (Delta Time)
        Output: Updates the Quatrenion vector
        '''
        Q_dot = 0.5 * (self._skewmatrix_wlb_b @ self._Quatrenion)
        self._Quatrenion =  self._Quatrenion + (delta_time * Q_dot)
        # print ("Updated Quatrenions",self._Quatrenion)
        self._Quatrenion = self._Quatrenion /  np.linalg.norm(self._Quatrenion)

    def UpdateRBL_Q (self):
        '''
        Input: None
        Ouput: Updates the Rotation from body to local level frame using the updated Quatrnions
        '''
         ##### First Row  #####
        self.rbl11 = ( self._Quatrenion[0]**2 ) - ( self._Quatrenion[1]**2 ) - ( self._Quatrenion[2]**2 ) + ( self._Quatrenion[3]**2 ) 

        self.rbl12 = 2 * ( ( self._Quatrenion[0] * self._Quatrenion[1] ) - ( self._Quatrenion[2] *  self._Quatrenion[3] ) )

        self.rbl13 = 2 * ( ( self._Quatrenion[0] * self._Quatrenion[2] ) + ( self._Quatrenion[1] * self._Quatrenion[3] ) ) 

        ##### Second Row #####
        self.rbl21 =  2 * ( ( self._Quatrenion[0] * self._Quatrenion[1] ) + ( self._Quatrenion[2] *  self._Quatrenion[3] ) )

        self.rbl22 =  - (self._Quatrenion[0]**2) + (self._Quatrenion[1]**2 ) - (self._Quatrenion[2]**2) + ( self._Quatrenion[3]**2 )

        self.rbl23 = 2 * ( ( self._Quatrenion[1] * self._Quatrenion[2] ) - ( self._Quatrenion[0] * self._Quatrenion[3] ) )

        ##### Third Row ######
        self.rbl31 = 2 * ( ( self._Quatrenion[0] * self._Quatrenion[2] ) - ( self._Quatrenion[1] * self._Quatrenion[3] ) )

        self.rbl32 = 2 * ( ( self._Quatrenion[1] * self._Quatrenion[2] ) + ( self._Quatrenion[0] * self._Quatrenion[3] ) )

        self.rbl33 = - (self._Quatrenion[0]**2) - (self._Quatrenion[1]**2 ) + (self._Quatrenion[2]**2) + ( self._Quatrenion[3]**2 )

        ##### Matrix formulation ######
        self._Rb_l = np.ndarray (shape= (3,3), dtype=float 
                    ,buffer= np.array ([[self.rbl11,self.rbl12,self.rbl13]
                    ,[self.rbl21,self.rbl22,self.rbl23]
                    ,[self.rbl31,self.rbl32,self.rbl33]]))
        # print ("Updated R_BL",self._Rb_l)

    def UpdateRBL_Attitude (self, roll, pitch, azimuth): 
        """
        Input: None, takes attitude components from the class initialization
        Output: Returns the initial rotation matrix from the body frame to the local level frame
        """
        # Rotation from body frame to local level frame
        ### Note: 1- The numpy trigonometric functions take radians as an input.
        ### 2- These matrix is intiated based on pre-saved entries / not real time. If real time the structure will have to change
        azimuth= np.deg2rad(azimuth)
        ##### First Row  #####
        self.rbl11 = ( np.cos(roll) * np.cos(azimuth) ) +  ( np.sin(roll) * np.sin(pitch) * np.sin(azimuth) ) 

        self.rbl12 = np.sin(azimuth) * np.cos(pitch)

        self.rbl13 = np.cos(azimuth)*np.sin(roll) - np.sin(azimuth)*np.sin(pitch)*np.cos(roll)
        ##### Second Row #####
        self.rbl21 = - np.sin(azimuth) * np.cos(roll) + np.cos(azimuth)*np.sin(pitch)*np.sin(roll)

        self.rbl22 = np.cos(azimuth)*np.cos(pitch) 

        self.rbl23 = - np.sin(azimuth)*np.sin(roll) - np.cos(azimuth)*np.sin(pitch)*np.cos(roll)
        ##### Third Row ######
        self.rbl31 = - np.cos(pitch)*np.sin(roll) 

        self.rbl32 = np.sin(pitch)

        self.rbl33 = np.cos(roll) * np.cos(pitch)
        ##### Matrix formulation ######
        self._Rb_l = np.ndarray (shape= (3,3), dtype=float 
                    ,buffer= np.array ([[self.rbl11,self.rbl12,self.rbl13]
                    ,[self.rbl21,self.rbl22,self.rbl23]
                    ,[self.rbl31,self.rbl32,self.rbl33]]))


    def UpdateAttitude (self):
        '''
        Input: None
        Output: Returns the updated attitude components (roll, pitch yaw (azimuth)) 
        '''
        self._pitch  = np.rad2deg(np.arctan2( self.rbl32 , np.sqrt( (self.rbl12)**2 + (self.rbl22)**2  ) ) )
        self._roll   = - np.rad2deg ( np.arctan2( self.rbl31 , self.rbl33 ) )
        self._azimuth    = np.rad2deg(np.arctan2( self.rbl12 , self.rbl22 ))
        # print ("attitude angles", self._roll, self._pitch ,self._azimuth)


    def OMEGA_IE_L (self):
        self._omega_ie_l = np.ndarray (shape= (3,3), dtype=float
                            ,buffer= np.array([[0, - self._wie_l[2], self._wie_l[1]]
                                ,[self._wie_l[2], 0 , - self._wie_l[0]]
                                ,[- self._wie_l[1], self._wie_l[0], 0]]) )
        # print ("OMEGA_IE_L: ",self._omega_ie_l)
    
    def OMEGA_EL_L (self):
        self._omega_el_l = np.ndarray (shape= (3,3), dtype=float 
                            ,buffer= np.array ([[0, - self._wel_l[2], self._wel_l[1]]
                                ,[self._wel_l[2], 0 , - self._wel_l[0]]
                                ,[- self._wel_l[1], self._wel_l[0], 0]]) )
        # print ("OMEGA_EL_L: ",self._omega_el_l)
    def UpdateAccelerometers (self, fx, fy, fz):
        self._fb = np.array([fx,fy,fz]).transpose()
        # print ("Accelerometers: ",self._fb)

    def UpdateG (self,g):
        self._gVector = np.array([ 0, 0, -g]).transpose() 
        # print ("updated g constant: ", self._gVector)
    def UpdateDeltaVelocity (self, delta_time):
        component_a = (self._Rb_l @ self._fb)
        # print ("delta velocity Component_a", component_a)
        component_b = ( (2 * self._omega_ie_l + self._omega_el_l) @ self._vl )
        # print ("delta velocity Component_a", component_b) 
        delta_v_t = component_a - component_b + (self._gVector)
        self._Delta_Vl_current = delta_v_t * delta_time
        # print ("self._Delta_Vl_current", self._Delta_Vl_current )
    
    def UpdateVelocity (self):
        self._vl = self._vl + 0.5 * (self._Delta_Vl_current + self._Delta_Vl_previous )
        # print ("Updated Velocities", self._vl )
        self._Delta_Vl_previous = self._Delta_Vl_current
        # self._vl = np.array([Ve,Vn,Vu])

    def UpdatePosition (self, prev_latitude, prev_longitude, prev_altitude, ve_Prev, vn_Prev, vu_Prev, delta_time):
        ve, vn, vu = self._vl[0], self._vl[1], self._vl[2]
        self._altitude = prev_altitude +  (0.5 * (vu + vu_Prev) * delta_time)
        self._latitude = prev_latitude + np.rad2deg ( (0.5* (vn + vn_Prev) * delta_time ) / (self._Rn + self._altitude) ) 
        self._longitude = prev_longitude + np.rad2deg ( ( 0.5 * (ve + ve_Prev)  * delta_time ) / ((self._Rn + prev_altitude) * np.cos(self._latitude)) ) 
        # print ("Updated Positions",  [self._latitude, self._longitude, self._altitude ])
         
    def CorrectAzimuth (self):
        if self._azimuth >= 360:
            self._azimuth = self._azimuth - 360.0
        elif self._azimuth < 0:
            self._azimuth = self._azimuth + 360.0
        # print ("Corrected azimuth", self._azimuth)

            
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
        self.WIE_L(latitude)
        self.WEL_L(latitude,altitude)
        self.WLB_B(wx,wy,wz)
        self.SkewMatrix_WLB_B()
        #### Stage 3: updating quatrenions ##################
        self.UpdateQuatrenion(dt)
        #### Stage 4: updating attitude components ########
        self.UpdateRBL_Q()
        self.UpdateAttitude()
        self.CorrectAzimuth ()
        ##### Stage 5: updating velocity components #######
        self.OMEGA_IE_L()
        self.OMEGA_EL_L()
        self.UpdateAccelerometers(fx,fy,fz)
        self.UpdateG(g)
        self.UpdateDeltaVelocity (dt)
        self.UpdateVelocity ()
        ##### Stage 6: updating position components #######
        self.UpdatePosition (latitude, longitude, altitude,  ve_Prev, vn_Prev, vu_Prev, dt)

        ##### Return the updated full navigation components #####
        return self._latitude, self._longitude, self._altitude , self._vl[0],self._vl[1], self._vl[2], self._roll, self._pitch, self._azimuth
                



       

