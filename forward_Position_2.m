% function to calculate the forward kinematics
%
%        X = forward_Position( q )
%
% where q is a column vector of joint angles

function X = forward_Position_2( q, a1, a2 )

% global a1 a2

X(1) = a1 * cos(q(1)) + a2 * cos(q(1) + q(2));

X(2) = a1 * sin(q(1)) + a2 * sin(q(1) + q(2));

X = X';
