% function to get angular displacement from position
%       q = inv_disp(X)
% where x is a column vector (2x1)

function q = inv_Position_2(X, a1, a2)

%global a1 a2

x = X(1);
y = X(2);

D = (x^2 + y^2 - a1^2 - a2^2)/(2*a1*a2);

q(2) = atan2(abs(sqrt(1-D^2)),D);
q(1) = atan2(y,x) - atan2(a2*sin(q(2)),(a1 + a2*cos(q(2))));

q = q';

