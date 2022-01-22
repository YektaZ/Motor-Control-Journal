% function to calculate the min jerk velocity at time t
%
%              Xd_dot = min_Jerk_Velocity( t, to, tf, Xo, Xf )
%
% where t is the current time, to is the start time for the
% min jerk trajectory and tf is the final time
% Xo is a column vector of the originating coordinates and
% Xf is a column vector of the final coordinates

function Xd_dot = min_Jerk_Velocity( t, to, tf, Xo, Xf )

T = (t-to)/(tf-to);
tfo = tf - to;

if (t <= to)
   Xd_dot = [0,0];
elseif ( (t > to) & (t <= tf) )
   Xd_dot(1) = (Xo(1) - Xf(1)) * ( 60*T^3 - 30*T^4 - 30*T^2 )/tfo;

   Xd_dot(2) = (Xo(2) - Xf(2)) * ( 60*T^3 - 30*T^4 - 30*T^2 )/tfo;
else
   Xd_dot = [0,0];
end

Xd_dot = Xd_dot';
