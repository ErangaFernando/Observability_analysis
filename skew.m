function ax = skew(a)
%% Produce the skew symmetric matrix from a 3x1 matrix
if ((size(a,1) == 3 || size(a,2) == 3) && numel(a) == 3)
    ax = [   0  -a(3)  a(2);
           a(3)    0  -a(1);
          -a(2)  a(1)    0];
end
    

end