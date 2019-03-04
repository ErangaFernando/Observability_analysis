% Housekeeping
clc;
close all;
clear all;

%% Variables
syms x y z qx qy qz qw bx by bz vx vy vz T real % states
syms wx wy wz ft real % Inputs
syms kpa kpp real % drag coefficients
syms m g t real 
syms p1x p1y p1z p2x p2y p2z p3x p3y p3z real
syms p4x p4y p4z real
syms d real
% Basis vector
e1 = [1 0 0]';
e2 = [0 1 0]';
e3 = [0 0 1]';



%% Trajectories & Scenarios
% traj 1.1
% 1m/s velocity in X direction
% Beacons perpendicular to the trajectory on the same plane

traj = 1.1;

if(traj == 1.1)
    q = [0.00000000000000e+000    5.04567886426745e-003    0.00000000000000e+000   999.987270481379e-003];
    q = [qx qy qz qw];
    % Position
    p = [x p1y p1z]';
    % Range beacon locations
    p1 = [p1x p1y p1z]';
    p2 = [0 p1y p1z]';
    p3 = [p1x p1y+d p1z]';
    % Velocity
    vb = [vx 0 0]'; % Body velocity
    K = diag([kpp,kpp,kpa]); % Drag matrix
    w = [0 0 0 0]'; % Gyroscope rates
end

out = observability_check( p,q,vb,w,p1,p2,p3,1)

