clear all;
clc;

syms x y z qx qy qz qw bx by bz vx vy vz T real % states
syms wx wy wz ft real % Inputs
syms kpa kpp real % drag coefficients
syms m g t real
syms p1x p1y p1z p2x p2y p2z p3x p3y p3z real
syms p4x p4y p4z real
syms ph th ps real
syms u real
syms k v real

% Basis vector
e1 = [1 0 0]';
e2 = [0 1 0]';
e3 = [0 0 1]';

p1 = [p1x p1y p1z]';
p2 = [p2x p2y p2z]';
p3 = [p3x p3y p3z]';
p = [x y z]';
eul = [ph th ps]';


K = diag([kpp,kpp,kpa]);
q = [qx qy qz qw]';
rpy = [ph th ps]';
Cq = quat2rot(q);
Rq = quatM_R([qx qy qz qw])/2;
Cr = rotsym('z',ps)*rotsym('y',th)*rotsym('x',ph);
Q =[[ qx,  qy,  qz,  qw];
    [ qy, -qx,  qw, -qz];
    [ qz, -qw, -qx,  qy]];
V =[[ vx,  vy,  vz,   0];
    [ vy, -vx,   0, -vz];
    [ vz,   0, -vx,  vy];
    [  0, -vz,  vy,  vx]];
G =[[ 0,  0,  g,   0];
    [ 0,  0,  0,  -g];
    [ g,  0,  0,   0];
    [ 0, -g,  0,   0]];
Qt = [[ qx,  -qy,  qz,  -qw];
      [ qy, qx, -qw,  -qz];
      [ qz,  -qw, -qx, qy]];

r1 = (p-p1);
r2 = (p-p2);
r3 = (p-p3);
Omega = [[1      tan(th)*sin(ph)     tan(th)*cos(ph)];
         [0              cos(ph)            -sin(ph)];
         [0      sin(ph)/cos(th)     cos(ph)/cos(th)]];
vb = [vx vy vz]';

X1 = [x y z vx vy vz qx qy qz qw T]';
X2 = [x y z vx vy vz ph th ps T]';

f0_1 = [Cq*vb; (-K*vb - T*e3 + Cq'*g*e3) ; zeros(4,1); 0];
f0_2 = [Cr*vb; (-K*vb - T*e3 + Cr'*g*e3) ; zeros(3,1); 0];
f1_1 = [zeros(3,4); zeros(3,4); Rq; zeros(1,4)];
f1_2 = [zeros(3,3); zeros(3,3); Omega; zeros(1,3)];
f2_1 = [zeros(3,1); zeros(3,1); zeros(4,1); 1];
f2_2 = [zeros(3,1); zeros(3,1); zeros(3,1); 1];

%% 3 range measurements
if(1)
    % Quaternion
    h = [-K*vb - T*e3; (p-p1)'*(p-p1)/2; (p-p2)'*(p-p2)/2; (p-p3)'*(p-p3)/2];
    for i = 1:size(h,1)
        for j = 1:size(X1,1)
            H1(i,j) = diff(h(i),X1(j));
        end
    end
    
    L1f0h = H1*f0_1;
    
    if(size(L1f0h,2)~=1)
        L1f0h = L1f0h(:);
    end
    for i = 1:size(L1f0h,1)
        for j = 1:size(X1,1)
            GL1f0h(i,j) = diff(L1f0h(i),X1(j));
        end
    end
    O1 = [H1;GL1f0h];
    
    %Gauss Elimination for three range measurements
    O46 = simplify(inv(O1(4:6,1:3))*O1(4:6,:));
    O1(4:6,:) = O46;
    O79 = inv(K)*O1(7:9,:) + O1(1:3,:);
    O1(7:9,:) = O79;
    O1012 = O1(10:12,:) - [vb vb vb]'*Cq'*O1(4:6,:);
    O1(10:12,:) = simplify(inv([r1 r2 r3]')*O1012);
    O1(1:3,:) = -inv(K)*O1(1:3,:);
    O1(10:12,:) = -Cq*O1(1:3,:) + O1(10:12,:)
end

%% Three range measurements with second order lie derivative
if(0)
    % Quaternion
    h = [-K*vb - T*e3; (p-p1)'*(p-p1)/2; (p-p2)'*(p-p2)/2; (p-p3)'*(p-p3)/2];
    for i = 1:size(h,1)
        for j = 1:size(X1,1)
            H1(i,j) = diff(h(i),X1(j));
        end
    end
    
    L1f0h = H1*f0_1;
    
    if(size(L1f0h,2)~=1)
        L1f0h = L1f0h(:);
    end
    for i = 1:size(L1f0h,1)
        for j = 1:size(X1,1)
            GL1f0h(i,j) = diff(L1f0h(i),X1(j));
        end
    end
    % Second order
    %L2f0f0h = GL1f0h*f0_1;
    % This give a mathematically true expression. But there is another
    % simpler way of representing that without most of the quaternions
    L2f0f0h = [K*K*f0_1(4:6,1); [[vb vb vb]'*vb + [r1 r2 r3]'*(Cq*(-K*vb - T*e3)  + g*e3 )]]
    for i = 1:size(L2f0f0h,1)
        for j = 1:size(X1,1)
            GL2f0f0h(i,j) = diff(L2f0f0h(i),X1(j));
        end
    end
    O1 = [H1;GL1f0h;GL2f0f0h];
    
    %Gauss Elimination for three range measurements
     O46 = simplify(inv(O1(4:6,1:3))*O1(4:6,:));
     O1(4:6,:) = O46;
     O1315 = O1(13:15,:) + K*O1(7:9,:);
     O1(13:15,:) = O1315;
     O79 = inv(K)*O1(7:9,:) + O1(1:3,:);
     O1(7:9,:) = O79;
     O1(1:3,:) = -inv(K)*O1(1:3,:);
     O1012 = O1(10:12,:) - [vb vb vb]'*Cq'*O1(4:6,:);
     O1(10:12,:) = simplify(inv([r1 r2 r3]')*O1012);
     O1(10:12,:) = -Cq*O1(1:3,:) + O1(10:12,:);
     O1618 = -[f0_1(4:6,:) f0_1(4:6,:) f0_1(4:6,:)]'*Cq'*O1(4:6,:) - (2*[vb vb vb]' - [r1 r2 r3]'*Cq*K)*O1(1:3,:) + O1(16:18,:);
     O1618(:,3) = zeros(3,1);
     O1(16:18,:) = O1618;
     simplify(inv([r1 r2 r3]')*O1618)
     simplify(O1)
end

%% Three range measurements and Euler angles
if(0)
    % Euler angles
    h = [-K*vb - T*e3; (p-p1)'*(p-p1)/2; (p-p2)'*(p-p2)/2; (p-p2)'*(p-p3)/2];
    for i = 1:size(h,1)
        for j = 1:size(X2,1)
            H12(i,j) = diff(h(i),X2(j));
        end
    end
    
    L1f0h2 = H12*f0_2;
    
    if(size(L1f0h2,2)~=1)
        L1f0h2 = L1f0h2(:);
    end
    for i = 1:size(L1f0h2,1)
        for j = 1:size(X2,1)
            GL1f0h2(i,j) = diff(L1f0h2(i),X2(j));
        end
    end
    O1 = [H12;GL1f0h2];
    % Euler angles
    O46 = simplify(inv(O1(4:6,1:3))*O1(4:6,:));
    O1(4:6,:) = O46;
    O79 = inv(K)*O1(7:9,:) + O1(1:3,:);
    O1(7:9,:) = O79;
    O1012 = O1(10:12,:) - [vb vb vb]'*Cr'*O1(4:6,:);
    O1(10:12,:) = simplify(inv([r1 r2 r3]')*O1012);
    O1(1:3,:) = -inv(K)*O1(1:3,:);
    O1(10:12,:) = -Cr*O1(1:3,:) + O1(10:12,:)
end
%% Two range measurements 
if(0)
    % Quaternion
    h = [-K*vb - T*e3; (p-p1)'*(p-p1)/2; (p-p2)'*(p-p2)/2];
    % Zeroth Order
    for i = 1:size(h,1)
        for j = 1:size(X1,1)
            H1(i,j) = diff(h(i),X1(j));
        end
    end
    % First order
    L1f0h = H1*f0_1;
    
    if(size(L1f0h,2)~=1)
        L1f0h = L1f0h(:);
    end
    for i = 1:size(L1f0h,1)
        for j = 1:size(X1,1)
            GL1f0h(i,j) = diff(L1f0h(i),X1(j));
        end
    end
    
    % Second order
    %L2f0f0h = GL1f0h*f0_1;
    % This give a mathematically true expression. But there is another
    % simpler way of representing that without most of the quaternions
    L2f0f0h = [K*K*f0_1(4:6,1); [[vb vb]'*vb + [r1 r2]'*(Cq*(-K*vb - T*e3)  + g*e3 )]];
    for i = 1:size(L2f0f0h,1)
        for j = 1:size(X1,1)
            GL2f0f0h(i,j) = diff(L2f0f0h(i),X1(j));
        end
    end
    
    O1 = [H1;GL1f0h;GL2f0f0h];
    
    %Gauss Elimination for three range measurements
    R2 = [r1 r2]';
    V2 = [vb vb]';
    VV = test2((-K*vb - T*e3));
    Ft = [f0_1(4:6) f0_1(4:6)]';
    O68 = inv(K)*O1(6:8,:) + O1(1:3,:);
    O1(6:8,:) = O68;
    O1113 = inv(K)*inv(K)*O1(11:13,:) - O1(1:3,:);
    O1(11:13,:) = O1113;
    O1(1:3,:) = -inv(K)*O1(1:3,:);
    O910 = O1(9:10,:) - R2*Cq*O1(1:3,:);
    O1(9:10,:) = O910;
    O1415 = O1(14:15,:) -(2*V2 - R2*Cq*K)*O1(1:3,:);
    O1(14:15,:) = O1415;
    O1113 = O1(11:13,:) + O1(6:8,:);
    O1(11:13,:) = O1113;
    O1(6:8,:) = -O1(6:8,:);
    O1 = simplify(O1);
    
    O_math = [zeros(3) eye(3) zeros(3,4) inv(K)*e3;
R2 zeros(2,3) zeros(2,4) zeros(2,1);
zeros(3) zeros(3) 2*Qt*G zeros(3,1);
V2*Cq' zeros(2,3) 2*R2*Q*V -R2*Cq*inv(K)*e3;
zeros(3,11);
Ft*Cq' zeros(2,3) 2*R2*Q*VV -2*V2*inv(K)*e3]
end


% for i = 1:size(h,1)
%     for j = 1:size(X2,1)
%         H2(i,j) = diff(h(i),X2(j));
%     end
% end
% 
% CqVb = Cq*vb;
% for i = 1:size(CqVb,1)
%     for j = 1:size(q,1)
%         dCqVb(i,j) = diff(CqVb(i),q(j));
%     end
% end
% 
% 
% r1CqVb = r1'*Cq*vb;
% for i = 1:size(r1CqVb,1)
%     for j = 1:size(q,1)
%         dr1CqVb(i,j) = diff(r1CqVb(i),q(j));
%     end
% end
% 
% Cqge3 = Cq'*g*e3;
% for i = 1:size(Cqge3,1)
%     for j = 1:size(q,1)
%         dCqge3(i,j) = diff(Cqge3(i),q(j));
%     end
% end
% 
% KCqge3 = K*Cq'*g*e3;
% for i = 1:size(KCqge3,1)
%     for j = 1:size(q,1)
%         dKCqge3(i,j) = diff(KCqge3(i),q(j));
%     end
% end
% 
% dKCqge3 - K*dCqge3
% for j = 1:size(q,1)
%     for i = 1:size(Cq,1)
%         for k = 1:size(Cq,2)          
%             dCqX(i,(j-1)*3+k) = diff(Cq(i,k),q(j));
%         end
%     end
% end
% 
% r1Cqvb = r1'*Cq*vb;
% for i = 1:size(r1Cqvb,1)
%     for j = 1:size(vb,1)
%         dr1Cqvb(i,j) = diff(r1Cqvb(i),vb(j));
%     end
% end
% 
% L1f0h = H1*f0_1;
% if(size(L1f0h,2)~=1)
%         L1f0h = L1f0h(:);
% end
% for i = 1:size(L1f0h,1)
%     for j = 1:size(X1,1)
%         GL1f0h(i,j) = diff(L1f0h(i),X1(j));
%     end
% end
% 
% L2f0f0h = GL1f0h*f0_1;
% 
% h22 = [-1/z 0 x/z^2; 0 -1/z y/z^2];
% h23 = h22(:);
% xb1 = [x y z]';
% for i = 1:size(h23,1)
%     for j = 1:size(xb1,1)
%         H23(i,j) = diff(h23(i),xb1(j));
%     end
% end
% 
% Crg = Cr*g*e3;
% 
% for i = 1:size(Crg,1)
%     for j = 1:size(eul,1)
%         dCrg(i,j) = diff(Crg(i),eul(j));
%     end
% end
% 
CrVb = Cr*vb;

for i = 1:size(CrVb,1)
    for j = 1:size(eul,1)
        dCrVb(i,j) = diff(CrVb(i),eul(j));
    end
end
% 
% Cqg = Cq'*g*e3;
% for i = 1:size(Cqg,1)
%     for j = 1:size(q,1)
%         dCqg(i,j) = diff(Cqg(i),q(j));
%     end
% end