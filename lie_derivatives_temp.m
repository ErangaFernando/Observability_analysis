%% This function will calcute the lie derivative

clear all;
clc;

%% model_no:
% 1: X = [x y z qx qy qz qw bx by bz vx vy vz T]';
% 2: X = [qx qy qz qw vx vy vz]';
% 3: X = [qx qy qz qw vx vy vz T]';
% 4: X = [x y z qx qy qz qw vx vy vz T]';
% 5: X = [x y z qx qy qz qw vx vy vz T]'; 4 positions available
model_no = 5;
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
if (model_no == 1)
    % State vector
    X = [x y z qx qy qz qw bx by bz vx vy vz T]';
    p = [x y z]'; % Position
    p1 = [p1x p1y p1z]';
    p2 = [p2x p2y p2z]';
    p3 = [p3x p3y p3z]';
    p4 = [p4x p4y p4z]';
    vb = [vx vy vz]'; % Body velocity
    Rq = quatM_R([qx qy qz qw]);
    bg = [bx by bz 0]'; %Gyroscope bias
    w = [wx wy wz 0]'; % Gyroscope rates
    K = diag([kpp,kpp,kpa]); % Drag matrix
    f_t = - K*vb/m + C'*e3*g;
    
    f0 = [C*vb; -Rq*bg/2; bg(1:3)/t; f_t - T*e3/m;0];
    f1 = [zeros(3,4); Rq; zeros(7,4)];
    f2 = [zeros(13,1); 1];
    
    
    h = [-K*vb/m - T*e3/m; e3'*p; (p'*p - p1'*p + p1'*p1 - p'*p1)/2];% (p'*p -p2'*p + p2'*p2 -p'*p2)/2];% (p'*p -p3'*p + p3'*p3 -p'*p3)/2];
    %H1 = zeros(size(h,1),size(X,1));
    for i = 1:size(h,1)
        for j = 1:size(X,1)
            H1(i,j) = diff(h(i),X(j));
        end
    end
    
    % 1st order lie derivatives
    L1f0h = H1*f0;
    L1f1h = H1*f1;
    L1f2h = H1*f2;
     
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
     
    % Second order lie derivaties
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
     
    GL2f0f0h
    GL2f0f1h
    GL2f0f2h
     
    O = [H1];
    O1 = [H1;GL1f0h];
    O2 = [H1;GL1f0h;GL2f0f0h];
    rank(O)
    rank(O1)
    rank(O2)

end

if (model_no == 2)
    % State vector
    X = [qx qy qz qw vx vy vz]';
    vb = [vx vy vz]'; % Body velocity
    Rq = quatM_R([qx qy qz qw]);
    w = [wx wy wz 0]'; % Gyroscope rates
    K = diag([kpp,kpp,kpa]); % Drag matrix
    f_t = - K*vb/m + C'*e3*g;
    
    f0 = [zeros(4,1); f_t];
    f1 = [Rq; zeros(3,4)];
    f2 = [zeros(6,1);1];
    
    h = -K*vb/m;
    
    % Zeroth Lie derivative
    % h
    for i = 1:size(h,1)
        for j = 1:size(X,1)
            H1(i,j) = diff(h(i),X(j));
        end
    end
    % Gradient of lie derivative : H1
    
    % 1st order lie derivatives
    L1f0h = H1*f0;
    L1f1h = H1*f1;
    L1f2h = H1*f2;
    
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
    
    % Second order lie derivaties
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
    
    GL2f0f0h
    GL2f0f1h
    GL2f0f2h
    
    O = [H1; GL1f0h; GL2f0f0h; GL2f0f1h; GL2f0f2h];    
    
end

if (model_no == 3)
    % State vector
    X = [qx qy qz qw vx vy vz T]';
    vb = [vx vy vz]'; % Body velocity
    Rq = quatM_R([qx qy qz qw]);
    w = [wx wy wz 0]'; % Gyroscope rates
    K = diag([kpp,kpp,kpa]); % Drag matrix
    f_t = - K*vb/m + C'*e3*g;
     
    f0 = [zeros(4,1); f_t - T*e3/m;0];
    f1 = [Rq; zeros(4,4)];
    f2 = [zeros(7,1); 1];
     
    h = -K*vb/m - T*e3/m;
     
    % Zeroth Lie derivative
    % h
    for i = 1:size(h,1)
        for j = 1:size(X,1)
            H1(i,j) = diff(h(i),X(j));
        end
    end
    % Gradient of lie derivative : H1
     
    % 1st order lie derivatives
    L1f0h = H1*f0;
    L1f1h = H1*f1;
    L1f2h = H1*f2;
     
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
     
    % Second order lie derivaties
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
     
    GL2f0f0h
    GL2f0f1h
    GL2f0f2h
     
    O = [H1; GL1f0h; GL2f0f0h; GL2f0f1h; GL2f0f2h];
     
end

if (model_no == 4)
    % State vector
    X = [x y z qx qy qz qw vx vy vz T]';
    p = [x y z]'; % Position
    p1 = [p1x p1y p1z]';
    p2 = [p2x p2y p2z]';
    p3 = [p3x p3y p3z]';
    p4 = [p4x p4y p4z]';
    vb = [vx vy vz]'; % Body velocity
    Rq = quatM_R([qx qy qz qw]);
    bg = [bx by bz 0]'; %Gyroscope bias
    w = [wx wy wz 0]'; % Gyroscope rates
    K = diag([kpp,kpp,kpa]); % Drag matrix
    f_t = - K*vb/m + C'*e3*g;
    
    f0 = [C*vb; zeros(4,1); f_t - T*e3/m;0];
    f1 = [zeros(3,4); Rq; zeros(4,4)];
    f2 = [zeros(10,1); 1];
    
    % Only height measurement
    h = [-K*vb/m - T*e3/m; e3'*p]; 
    %H1 = zeros(size(h,1),size(X,1));
    for i = 1:size(h,1)
        for j = 1:size(X,1)
            H1(i,j) = diff(h(i),X(j));
        end
    end
    
    % 1st order lie derivatives
    L1f0h = H1*f0;
    L1f1h = H1*f1;
    L1f2h = H1*f2;
     
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
     
    % Second order lie derivaties
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
     
    GL2f0f0h
    GL2f0f1h
    GL2f0f2h
     
    O = [H1; GL1f0h; GL2f0f0h; GL2f0f1h; GL2f0f2h];
     

end

if (model_no == 5)
    % State vector
    X = [x y z qx qy qz qw vx vy vz T]';
    p = [x y z]'; % Position
    p1 = [p1x p1y p1z]';
    p2 = [p2x p2y p2z]';
    p3 = [p3x p3y p3z]';
    p4 = [p4x p4y p4z]';
    vb = [vx vy vz]'; % Body velocity
    Rq = quatM_R([qx qy qz qw]);
    bg = [bx by bz 0]'; %Gyroscope bias
    w = [wx wy wz 0]'; % Gyroscope rates
    K = diag([kpp,kpp,kpa]); % Drag matrix
    f_t = - K*vb/m + C'*e3*g;
    
    f0 = [C*vb; zeros(4,1); f_t - T*e3/m;0];
    f1 = [zeros(3,4); Rq/2; zeros(4,4)];
    f2 = [zeros(10,1); 1];
    
    
    h = [-K*vb/m - T*e3/m; e3'*p; (p'*p - p1'*p + p1'*p1 - p'*p1)/2; (p'*p -p2'*p + p2'*p2 -p'*p2)/2];% (p'*p -p3'*p + p3'*p3 -p'*p3)/2];% (p'*p -p4'*p + p4'*p4 -p'*p4)/2];
    %H1 = zeros(size(h,1),size(X,1));
    for i = 1:size(h,1)
        for j = 1:size(X,1)
            H1(i,j) = diff(h(i),X(j));
        end
    end
    
    % 1st order lie derivatives
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
     
    % Second order lie derivaties
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
     
    O = [H1];
    O1 = [H1;GL1f0h];
    O2 = [H1;GL1f0h;GL2f0f0h];
    rank(O)
    rank(O1)
    rank(O2)
     
end