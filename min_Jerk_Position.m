% function to calculate the min jerk position at time t
%
%                 Xd = min_Jerk_Position( t, to, tf, Xo, Xf )
%
% where t is the current time, to is the start time for the
% min jerk trajectory and tf is the final time
% Xo is a column vector of the originating coordinates and
% Xf is a column vector of the final coordinates

function Xd = min_Jerk_Position( t, to, tf, Xo, Xf )

T = (t-to)/(tf-to);

if (t <= to)
   Xd = Xo';
elseif ( (t > to) & (t <= tf) )
   Xd(1) = Xo(1) + (Xo(1) - Xf(1)) * ( 15*T^4 - 6*T^5 - 10*T^3 );

   Xd(2) = Xo(2) + (Xo(2) - Xf(2)) * ( 15*T^4 - 6*T^5 - 10*T^3 );
else
   Xd = Xf';
end

Xd = Xd';
