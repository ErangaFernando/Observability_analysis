function V = test2(v)
vx = v(1); vy = v(2); vz = v(3);
V = [[ vx,  vy,  vz,   0];
[ vy, -vx,   0, -vz];
[ vz,   0, -vx,  vy];
[  0, -vz,  vy,  vx]]
end
