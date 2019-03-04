%% This function will calcute the lie derivatives
% Velocity in the global frame is estimated
% Z direction is up
%


clear all;
clc;

%% Model No 1
% Model 1: This seems to go nowhere... Have to move back to the velocity in the
% quadrotor frame 
%% Model No 2
% In this model the accelerometer and the thrust are taken as an input of the system. The
% measurement is the position of the quadrotor.
model_no = 2;


syms x y z qx qy qz qw bx by bz vx vy vz T real % states
syms wx wy wz ft real % Inputs
syms kpa kpp real % drag coefficients
syms m g t real 
syms p1x p1y p1z p2x p2y p2z p3x p3y p3z real
syms p4x p4y p4z real

% Basis vector
e1 = [1 0 0]';
e2 = [0 1 0]';
e3 = [0 0 1]';

% Rotation matrix
C = quat2rot([qx qy qz qw]); % Rotation matrix
%% Different models
if(model_no == 1)
    X = [x y z qx qy qz qw vx vy vz T]';
    p = [x y z]'; % Position
    p1 = [p1x p1y p1z]';
    p2 = [p2x p2y p2z]';
    p3 = [p3x p3y p3z]';
    p4 = [p4x p4y p4z]';
    v = [vx vy vz]'; % Body velocity
    Rq = quatM_R([qx qy qz qw]);
    bg = [bx by bz 0]'; %Gyroscope bias
    w = [wx wy wz 0]'; % Gyroscope rates
    K = diag([kpp,kpp,kpa]); % Drag matrix
    f_t = -C*K*C'*v/m - e3*g;
    
    %% Process model Affine form
    f0 = [v; zeros(4,1); f_t + T*C*e3/m;0];
    f1 = [zeros(3,4); Rq/2; zeros(4,4)];
    f2 = [zeros(10,1); 1];
    
    %% Measurement model
%     % Zero range measurements
%     h = [-K*v/m - T*e3/m; e3'*p];
    % One range measurement
    h = [-K*C'*v/m + T*e3/m; e3'*p; (p'*p - p1'*p + p1'*p1 - p'*p1)/2];
%     % Two range measurements
%     h = [-K*v/m - T*e3/m; e3'*p; (p'*p - p1'*p + p1'*p1 - p'*p1)/2 (p'*p -p2'*p + p2'*p2 -p'*p2)/2];
%     % Three range measurements
%     h = [-K*v/m - T*e3/m; e3'*p; (p'*p - p1'*p + p1'*p1 - p'*p1)/2 (p'*p -p2'*p + p2'*p2 -p'*p2)/2 (p'*p -p3'*p + p3'*p3 -p'*p3)/2];
%     
    for i = 1:size(h,1)
        for j = 1:size(X,1)
            H1(i,j) = diff(h(i),X(j));
        end
    end
    
    %% 1st order lie derivatives
    L1f0h = H1*f0;
    L1f1h = H1*f1;
    L1f2h = H1*f2;
    
    size(L1f0h)
     
    if(size(L1f0h,2)~=1)
        L1f0h = L1f0h(:);
    end
    for i = 1:size(L1f0h,1)
        for j = 1:size(X,1)
            GL1f0h(i,j) = diff(L1f0h(i),X(j));
        end
    end
     
    if(size(L1f1h,2)~=1)
        L1f1h = L1f1h(:);
    end
    for i = 1:size(L1f1h,1)
        for j = 1:size(X,1)
            GL1f1h(i,j) = diff(L1f1h(i),X(j));
        end
    end
     
    if(size(L1f2h,2)~=1)
        L1f2h = L1f2h(:);
    end
    for i = 1:size(L1f2h,1)
        for j = 1:size(X,1)
            GL1f2h(i,j) = diff(L1f2h(i),X(j));
        end
    end
     
    %% Second order lie derivaties
    L2f0f0h = GL1f0h*f0;
    L2f0f1h = GL1f0h*f1;
    L2f0f2h = GL1f0h*f2;
     
    if(size(L2f0f0h,2)~=1)
        L2f0f0h = L2f0f0h(:);
    end
    for i = 1:size(L2f0f0h,1)
        for j = 1:size(X,1)
            GL2f0f0h(i,j) = diff(L2f0f0h(i),X(j));
        end
    end
     
     
    if(size(L2f0f1h,2)~=1)
        L2f0f1h = L2f0f1h(:);
    end
    for i = 1:size(L2f0f1h,1)
        for j = 1:size(X,1)
            GL2f0f1h(i,j) = diff(L2f0f1h(i),X(j));
        end
    end
     
     
    if(size(L2f0f2h,2)~=1)
        L2f0f2h = L2f0f2h(:);
    end
    for i = 1:size(L2f0f2h,1)
        for j = 1:size(X,1)
            GL2f0f2h(i,j) = diff(L2f0f2h(i),X(j));
        end
    end
    
    
    %% Third order lie derivatives
    L3f0f0f0h = GL2f0f0h*f0;
    L3f0f0f1h = GL2f0f0h*f1;
    L3f0f0f2h = GL2f0f0h*f2; %Constant
    
    L3f0f1f0h = GL2f0f1h*f0;
    L3f0f1f1h = GL2f0f1h*f1; 
    L3f0f1f2h = GL2f0f1h*f2; % Only zero
    
    
    if(size(L3f0f0f0h,2)~=1)
        L3f0f0f0h = L3f0f0f0h(:);
    end
    for i = 1:size(L3f0f0f0h,1)
        for j = 1:size(X,1)
            GL3f0f0f0h(i,j) = diff(L3f0f0f0h(i),X(j));
        end
    end
    
    if(size(L3f0f0f1h,2)~=1)
        L3f0f0f1h = L3f0f0f1h(:);
    end
    for i = 1:size(L3f0f0f1h,1)
        for j = 1:size(X,1)
            GL3f0f0f1h(i,j) = diff(L3f0f0f1h(i),X(j));
        end
    end
    
    if(size(L3f0f0f2h,2)~=1)
        L3f0f0f2h = L3f0f0f2h(:);
    end
    for i = 1:size(L3f0f0f2h,1)
        for j = 1:size(X,1)
            GL3f0f0f2h(i,j) = diff(L3f0f0f2h(i),X(j));
        end
    end
    
    if(size(L3f0f1f0h,2)~=1)
        L3f0f1f0h = L3f0f1f0h(:);
    end
    for i = 1:size(L3f0f1f0h,1)
        for j = 1:size(X,1)
            GL3f0f1f0h(i,j) = diff(L3f0f1f0h(i),X(j));
        end
    end
    
    if(size(L3f0f1f1h,2)~=1)
        L3f0f1f1h = L3f0f1f1h(:);
    end
    for i = 1:size(L3f0f1f1h,1)
        for j = 1:size(X,1)
            GL3f0f1f1h(i,j) = diff(L3f0f1f1h(i),X(j));
        end
    end
    
    if(size(L3f0f1f2h,2)~=1)
        L3f0f1f2h = L3f0f1f2h(:);
    end
    for i = 1:size(L3f0f1f2h,1)
        for j = 1:size(X,1)
            GL3f0f1f2h(i,j) = diff(L3f0f1f2h(i),X(j));
        end
    end
    
    %%
    O = [H1];
    O1 = [H1;GL1f0h];
    O2 = [H1;GL1f0h;GL2f0f0h];
    O3 = [H1;GL1f0h;GL2f0f1h];
    O4_1 = [H1;GL1f0h;GL2f0f1h;GL3f0f0f0h];
    O4_2 = [H1;GL1f0h;GL2f0f1h;GL3f0f0f1h];
    O4_3 = [H1;GL1f0h;GL2f0f1h;GL3f0f0f2h];
    O4_4 = [H1;GL1f0h;GL2f0f1h;GL3f0f1f0h];
    O4_5 = [H1;GL1f0h;GL2f0f1h;GL3f0f1f1h];
    O4_6 = [H1;GL1f0h;GL2f0f1h;GL3f0f1f2h];
    O5_1 = [H1;GL1f0h;GL2f0f0h;GL3f0f0f0h];
    O5_2 = [H1;GL1f0h;GL2f0f0h;GL3f0f0f1h];
    O5_3 = [H1;GL1f0h;GL2f0f0h;GL3f0f0f2h];
    O5_4 = [H1;GL1f0h;GL2f0f0h;GL3f0f1f0h];
    O5_5 = [H1;GL1f0h;GL2f0f0h;GL3f0f1f1h];
    O5_6 = [H1;GL1f0h;GL2f0f0h;GL3f0f1f2h];
    r_O = rank(O)
    r_O1 = rank(O1)
    r_O2 = rank(O2)
    r_O3 = rank(O3)
    r_O4_1 = rank(O4_1)
    r_O4_2 = rank(O4_2)
    r_O4_3 = rank(O4_3)
    r_O4_4 = rank(O4_4)
    r_O4_5 = rank(O4_5)
    r_O4_6 = rank(O4_6)
    r_O5_1 = rank(O5_1)
    r_O5_2 = rank(O5_2)
    r_O5_3 = rank(O5_3)
    r_O5_4 = rank(O5_4)
    r_O5_5 = rank(O5_5)
    r_O5_6 = rank(O5_6)
end

if(model_no == 2)
    X = [x y z qx qy qz qw vx vy vz]';
    p = [x y z]'; % Position
    p1 = [p1x p1y p1z]';
    p2 = [p2x p2y p2z]';
    p3 = [p3x p3y p3z]';
    p4 = [p4x p4y p4z]';
    v = [vx vy vz]'; % Body velocity
    Rq = quatM_R([qx qy qz qw]);
    bg = [bx by bz 0]'; %Gyroscope bias
    w = [wx wy wz 0]'; % Gyroscope rates
    K = diag([kpp,kpp,kpa]); % Drag matrix
    f_t = -C*K*C'*v/m - e3*g;
    
    %% Process model Affine form
    Ki = inv(-K);
    f0 = [zeros(3,1); zeros(4,1); f_t];
    f1 = [C*Ki; zeros(7,3)];
    f2 = [C*Ki*e3; zeros(4,1); C*e3];
    f3 = [zeros(3,4); Rq/2; zeros(3,4)];
    
    %% Measurement model;
    h = [p; e3'*p];
    
    %% Zeroth order Lie derivatives
    for i = 1:size(h,1)
        for j = 1:size(X,1)
            H1(i,j) = diff(h(i),X(j));
        end
    end
    %% 1st order lie derivatives
    L1f0h = H1*f0; % All zeros
    L1f1h = H1*f1;
    L1f2h = H1*f2; 
    L1f3h = H1*f3; % All zeros
    
    if(size(L1f1h,2)~=1)
        L1f1h = L1f1h(:);
    end
    for i = 1:size(L1f1h,1)
        for j = 1:size(X,1)
            GL1f1h(i,j) = diff(L1f1h(i),X(j));
        end
    end
     
    if(size(L1f2h,2)~=1)
        L1f2h = L1f2h(:);
    end
    for i = 1:size(L1f2h,1)
        for j = 1:size(X,1)
            GL1f2h(i,j) = diff(L1f2h(i),X(j));
        end
    end
    
    %% Second order Lie derivatives
    L2f1f0h = GL1f1h*f0; % All zeros
    L2f1f1h = GL1f1h*f1; % All zeros
    L2f1f2h = GL1f1h*f2; % All zeros
    L2f1f3h = GL1f1h*f3; 
    
    L2f2f0h = GL1f2h*f0; % All zeros
    L2f2f1h = GL1f2h*f1; % All zeros
    L2f2f2h = GL1f2h*f2; % All zeros
    L2f2f3h = GL1f2h*f3;
    
    if(size(L2f1f3h,2)~=1)
        L2f1f3h = L2f1f3h(:);
    end
    if(size(L2f2f3h,2)~=1)
        L2f2f3h = L2f2f3h(:);
    end
    
    for i = 1:size(L2f1f3h,1)
        for j = 1:size(X,1)
            GL2f1f3h(i,j) = diff(L2f1f3h(i),X(j));
        end
    end
    
    for i = 1:size(L2f2f3h,1)
        for j = 1:size(X,1)
            GL2f2f3h(i,j) = diff(L2f2f3h(i),X(j));
        end
    end
    O1 = [H1];
    O2 = [H1; GL1f1h];
    O3 = [H1; GL1f2h];
    O4 = [H1; GL1f1h;GL2f1f3h];
    O5 = [H1; GL1f2h;GL2f1f3h];
    O6 = [H1; GL1f1h;GL2f2f3h];
    O7 = [H1; GL1f2h;GL2f2f3h];
    
end