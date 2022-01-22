% function to get angular accelerations from translational
% velocity and position,
% 	dq/dt = inv_Velocity(dX/dt, q)
% where X, theta and theta' are col vectors (2x1)

function q_dot = inv_Velocity(X_dot, q, a1, a2)

%global a1 a2

q1 = q(1);
q2 = q(2); 

J11 = a2 * cos( q1+q2 );
J12 = a2 * sin( q1+q2 );
J21 = -( a1 * cos(q1) + a2 * cos( q1+q2 ) );
J22 = -( a1 * sin(q1) + a2 * sin( q1+q2 ) );

J_inv = [J11 J12;
	J21 J22]/( a1 * a2 * sin(q2) );

q_dot = J_inv * X_dot;

