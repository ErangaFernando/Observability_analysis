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


f0 = [Cq*vb; (-K*vb - T*e3 + Cq'*g*e3) ; zeros(4,1); 0];
f1 = [zeros(3,4); zeros(3,4); Rq; zeros(1,4)];
f2 = [zeros(3,1); zeros(3,1); zeros(4,1); 1];


%% 3 range measurements

% Quaternion
h = [-K*vb - T*e3; (p-p1)'*(p-p1)/2; (p-p2)'*(p-p2)/2; (p-p3)'*(p-p3)/2];
for i = 1:size(h,1)
    for j = 1:size(X1,1)
        H1(i,j) = diff(h(i),X1(j));
    end
end

L1f0h = H1*f0;

if(size(L1f0h,2)~=1)
    L1f0h = L1f0h(:);
end
for i = 1:size(L1f0h,1)
    for j = 1:size(X1,1)
        GL1f0h(i,j) = diff(L1f0h(i),X1(j));
    end
end
O1 = [H1;GL1f0h];

b1 = h;

for i = 1:size(b1,1)
    for j = 1:size(X1,1)
        db1(i,j) = diff(b1(i),X1(j));
    end
end

db1f0 = db1*f0;
db1f1 = db1*f1;
db1f2 = db1*f2;

b11 = [-K*e1; zeros(3,1)];
b12 = [-K*e2; zeros(3,1)];
b13 = [-K*e3; zeros(3,1)];

b2 = [-K*Cq'*g*e3; zeros(3,1)];
b3 = [zeros(3,1); [r1 r2 r3]'*Cq*vb];

for i = 1:size(b11,1)
    for j = 1:size(X1,1)
        db11(i,j) = diff(b11(i),X1(j));
    end
end

for i = 1:size(b12,1)
    for j = 1:size(X1,1)
        db12(i,j) = diff(b12(i),X1(j));
    end
end

for i = 1:size(b13,1)
    for j = 1:size(X1,1)
        db13(i,j) = diff(b13(i),X1(j));
    end
end

for i = 1:size(b2,1)
    for j = 1:size(X1,1)
        db2(i,j) = diff(b2(i),X1(j));
    end
end

for i = 1:size(b3,1)
    for j = 1:size(X1,1)
        db3(i,j) = diff(b3(i),X1(j));
    end
end

db3*f0