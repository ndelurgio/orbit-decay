%--------------------------------------------------------------------------
%
% Ephemeris computation using variable-order Radau IIA integrator with
% step-size control
%
% Last modified:   2018/02/11   Meysam Mahooti
%
%--------------------------------------------------------------------------
function Eph = Ephemeris(Y0, N_Step, Step)

options = rdpset('RelTol',1e-8,'AbsTol',1e-10);
% [t,yout] = radau(@Accel,(0:Step:N_Step*Step),Y0,options);
[t,yout] = radau(@Accel,(0:Step:N_Step*Step),Y0,options);
Eph(:,1) = t;
Eph(:,2:7) = yout;

