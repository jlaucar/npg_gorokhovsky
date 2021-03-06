function[x_new, time_new] = lorenz_96_adv_1step(x, time, forcing)
%% lorenz_96_adv_1step Does a single time step advance for lorenz_96 40-variable model using four step runge-kutta time step
%
% x is the 40-vector state, time is the 2-vector days and seconds time

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id: lorenz_96_adv_1step.m 11694 2017-06-02 21:39:25Z nancy@ucar.edu $

global DELTA_T

% Add in ability to return the raw dt for getting the slope if time is sent as negative
if(time < 0) 
   x_new = comp_dt(x, forcing);
   time_new = time + 1;
else

   % Compute first intermediate step
   dx = comp_dt(x, forcing);
   x1 = DELTA_T * dx;
   inter = x + x1 / 2;

   % Compute second intermediate step
   dx = comp_dt(inter, forcing);
   x2 = DELTA_T * dx;
   inter = x + x2 / 2;

   % Compute third intermediate step
   dx = comp_dt(inter, forcing);
   x3 = DELTA_T * dx;
   inter = x + x3;

   % Compute fourth intermediate step
   dx = comp_dt(inter, forcing);
   x4 = DELTA_T * dx;

   % Compute new value for x
   x_new = x + x1/6 + x2/3 + x3/3 + x4/6;

   % Increment time step
   time_new = time + 1;
end

end

%------------------------------------------------------------------------------

function[dt] = comp_dt(x,forcing)

global MODEL_SIZE

dt = zeros(1,MODEL_SIZE);

for j = 1:MODEL_SIZE
   jp1 = j + 1;
   if(jp1 > MODEL_SIZE), jp1 = 1; end

   jm2 = j - 2;
   if(jm2 < 1), jm2 = MODEL_SIZE + jm2; end

   jm1 = j - 1;
   if(jm1 < 1), jm1 = MODEL_SIZE; end

   dt(j) = (x(jp1) - x(jm2)) * x(jm1) - x(j) + forcing;
end

end

% <next few lines under version control, do not edit>
% $URL: https://svn-dares-dart.cgd.ucar.edu/DART/releases/Manhattan/documentation/DART_LAB/matlab/private/lorenz_96_adv_1step.m $
% $Revision: 11694 $
% $Date: 2017-06-02 15:39:25 -0600 (Fri, 02 Jun 2017) $
