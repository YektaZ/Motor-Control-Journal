% function to calculate the min jerk acceleration at time t
%
%          Xd_ddot = min_Jerk_Acceleration( t, to, tf, Xo, Xf )
%
% where t is the current time, and tf is the final time
% Xo is a column vector of the originating coordinates and
% Xf is a column vector of the final coordinates

function Xd_ddot = min_Jerk_Acceleration( t, to, tf, Xo, Xf )

T = (t-to)/(tf-to);
tfo = tf - to;

if (t <= to)
   Xd_ddot = [0,0];
elseif ( (t > to) & (t <= tf) )
   Xd_ddot(1) = (Xo(1) - Xf(1)) * ( 180*T^2 - 120*T^3 - 60*T )/(tfo*tfo);

   Xd_ddot(2) = (Xo(2) - Xf(2)) * ( 180*T^2 - 120*T^3 - 60*T )/(tfo*tfo);
else
   Xd_ddot = [0,0];
end

Xd_ddot = Xd_ddot';
