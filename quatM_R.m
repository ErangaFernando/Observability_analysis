function R = quatM_R(q)
%% This function will produce the right multiplication matrix for the quaternion q
% Input : q : A Quaternion, can be either a column vector or a row vector
% Output: L : Right multiplication matrix of q

if(numel(q) ~= 4)
    error('Check dimensions of the quaternion');
end

if (size(q,2) == 4)
    q = q';
end

q_x = [   0  -q(3)  q(2);
        q(3)    0  -q(1);
       -q(2)  q(1)    0 ];
   
R = [(q(4)*eye(3) + q_x) q(1:3,1);
          -q(1:3,1)'      q(4)];

end