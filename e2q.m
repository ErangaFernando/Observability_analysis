function q = e2q(roll,pitch,yaw)
%% This function will genereate the quaternion q related to R-P-Y Euler angles
% Output: q : [q_x q_y q_z q_w]'
% Input : roll, pitch, yaw : Angles are in radians and rotation is Z-Y-X
%       
Q23 = [sin(roll/2) 0 0 cos(roll/2)]';
Q12 = [0 sin(pitch/2) 0 cos(pitch/2)]';
Qw1 = [0 0 sin(yaw/2) cos(yaw/2)]';

q = quatM_L(Q23)*quatM_L(Q12)*Qw1;

end