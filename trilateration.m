clear all;
clc

for count = 1: 10000

point = 10*rand(1)*rand(3,1);
p11 = 10*rand(1)*rand(3,1);
p21 = 10*rand(1)*rand(3,1);
p31 = 10*rand(1)*rand(3,1);



r1 = norm(point-p11);
r2 = norm(point-p21);
r3 = norm(point-p31);

%Calculating the transformation 
zb = cross((p21-p11),(p31-p11));
zb = zb/norm(zb);
xb = (p21-p11)/norm(p21-p11);
yb = cross(zb,xb);

R = [xb yb zb];
T = [R p11; 0 0 0 1];
Tinv  = inv(T);

p1n = Tinv*[p11;1];
p2n = Tinv*[p21;1];
p3n = Tinv*[p31;1];
pointn = Tinv*[point;1];

r1 - norm(p1n-pointn);
r2 - norm(p2n-pointn);
r3 - norm(p3n-pointn);

p1 = p1n(1:3);
p2 = p2n(1:3);
p3 = p3n(1:3);

p1'
p = p2'
p3'







%% Calculating the coordinate frame
ex = (p2-p1)/norm(p2-p1);
i = dot(ex,(p3-p1));
ey = (p3 - p1 - i*ex)/norm(p3 - p1 - i*ex);
ez = cross(ex,ey);

d = norm(p2-p1);

j = dot(ey,(p3 - p1));

%% Calculating the xyz
x = (r1^2 - r2^2 + d^2)/(2*d);
y = (r1^2 - r3^2 + i^2 + j^2)/(2*j) - i*x/j;
z = sqrt(r1^2 - x^2 - y^2);

point_cal1 = x*ex + y*ey + z*ez;
point_cal2 = x*ex + y*ey - z*ez;

point_cal1 = T*[point_cal1;1];
point_cal1 = point_cal1(1:3);
point_cal2 = T*[point_cal2;1];
point_cal2 = point_cal2(1:3);


err = (point_cal1 - point)'
err2 = (point_cal2 - point)'

% if(abs(err(1))>1.0e-10 || abs(err(2))>1.0e-10 || abs(err(3))>1.0e-10)
%     [p1(3) p2(3) p3(3)];
%     rank([p11 p21 p31])
%     count
%     err
%     
%     d < r1+r2
%     d> r2-r1
%     error('Error is too large')
% end



end