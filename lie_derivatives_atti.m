clear all;
clc;

%% model_no:
% 1: X = [x y z qx qy qz qw bx by bz vx vy vz T]';

model_no = 5;
syms x y z qx qy qz qw bx by bz vx vy vz T real % states
syms wx wy wz ft real % Inputs
syms kpa kpp real % drag coefficients
syms m g t real 
syms p1x p1y p1z p2x p2y p2z p3x p3y p3z real
e3 = [0 0 1]';

% Rotation matrix
C = quat2rot([qx qy qz qw]); % Rotation matrix

X = [qx qy qz qw]';
Rq = quatM_R([qx qy qz qw]);
w = [wx wy wz 0]'; % Gyroscope rates

f1 = Rq*w/2;
h = -C'*g*e3;

for i = 1:size(h,1)
        for j = 1:size(X,1)
            H1(i,j) = diff(h(i),X(j));
        end
end

 % 1st order lie derivatives
    L1f1h = H1*f1;
    
    for i = 1:size(L1f1h,1)
        for j = 1:size(X,1)
            GL1f1h(i,j) = diff(L1f1h(i),X(j));
        end
    end
    
    