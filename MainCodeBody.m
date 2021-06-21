% Initialization
Rb_l = [   cos(roll(1))*cos(yaw(1))-sin(roll(1))*sin(pitch(1))*sin(yaw(1))       -cos(pitch(1))*sin(yaw(1))	  sin(roll(1))*cos(yaw(1))+cos(roll(1))*sin(pitch(1))*sin(yaw(1)) ;
         cos(roll(1))*sin(yaw(1))+sin(roll(1))*sin(pitch(1))*cos(yaw(1))    cos(pitch(1))*cos(yaw(1))   sin(roll(1))*sin(yaw(1))-cos(roll(1))*sin(pitch(1))*cos(yaw(1)) ;
                    -sin(roll(1))*cos(pitch(1))                             sin(pitch(1))                     cos(roll(1))*cos(pitch(1))                      ];

q4 = sqrt(1+Rb_l(1,1)+Rb_l(2,2)+Rb_l(3,3))/2;
q1= 0.25*(Rb_l(3,2)-Rb_l(2,3))/q4;
q2 = 0.25*(Rb_l(1,3)-Rb_l(3,1))/q4;
q3 = 0.25*(Rb_l(2,1)-Rb_l(1,2))/q4;

Quaternion = [q1 q2 q3 q4]';
Quaternion = Quaternion/norm(Quaternion);

delta_vl_prev=[0 0 0]';

Local_g = a1*(1+(a2*sin(pi/180*BP_Lat(1))^2)+(a3*sin(pi/180*BP_Lat(1))^4))+(a4+a5*sin(pi/180*BP_Lat(1))^2)*BP_Alt(1)+a6*BP_Alt(1)^2;

% Loop i=1
% Function call

[latitude(i+1), longitude(i+1),altitude(i+1),ve(i+1),vn(i+1),vu(i+1), ...
            pitch(i+1), roll(i+1), azi(i+1),yaw(i+1), Rb_l,Quaternion, delta_vl_prev] ...
            = mech_FulI_IMU(sample_time(i), w.x(i), w.y(i), w.z(i), f.x(i),f.y(i),f.z(i),...
            latitude(i), longitude(i),altitude(i),ve(i),vn(i),vu(i), Rb_l, Quaternion, delta_vl_prev, Local_g);
% end Loop i
        