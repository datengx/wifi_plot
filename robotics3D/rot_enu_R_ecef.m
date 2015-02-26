function M = rot_enu_R_ecef(phi0, lambda0)
% Construct the matrix that rotates Cartesian vectors from geocentric to
% local vertical.

% Copyright 2005 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2005/11/15 01:37:07 $

sinphi0 = sin(phi0);
cosphi0 = cos(phi0);
sinlambda0 = sin(lambda0);
coslambda0 = cos(lambda0);

M = [0 1 0; 1 0 0; 0 0 -1] * [    1      0        0   ; ...
         0   sinphi0  cosphi0; ...
         0  -cosphi0  sinphi0] ...
  * [-sinlambda0   coslambda0  0; ...
     -coslambda0  -sinlambda0  0; ...
           0          0        1];
       
       
%M = rotz(lambda0) * rotx(phi0);