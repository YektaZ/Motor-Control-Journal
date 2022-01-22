% function to calculate forward velocity kinematics
%
%        X_dot = forward_Velocity( q, q_dot )
%
% where q and q_dot are column vectors of joint coordinates

function X_dot = forward_Velocity( q, q_dot, a1, a2 )

%global a1 a2

q1 = q(1);
q2 = q(2);

J11 = -a1 * sin(q1) - a2 * sin(q1+q2);
J12 = -a2 * sin(q1+q2);
J21 = a1 * cos(q1) + a2 * cos(q1+q2);
J22 = a2 * cos(q1+q2);

% X_dot = J*q_dot

X_dot(1) = J11 * q_dot(1) + J12 * q_dot(2);
X_dot(2) = J21 * q_dot(1) + J22 * q_dot(2);

X_dot = X_dot';
