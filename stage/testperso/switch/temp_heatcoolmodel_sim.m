function [xn, d, z, y] = temp_heatcoolmodel_sim(x, u, params)
% [xn, d, z, y] = temp_heatcoolmodel_sim(x, u, params)
% simulates the hybrid system one step ahead.
% Parameters:
%   x: current state
%   u: input
%   params: structure containing values for
%           all symbolic parameters
% Output:
%   xn: state in the next timestep
%   u: output
%   d, z: Boolean and real auxiliary variables
%
% HYSDEL 2.0.5 (Build: 20111112)
% Copyright (C) 1999-2002  Fabio D. Torrisi
% 
% HYSDEL comes with ABSOLUTELY NO WARRANTY;
% HYSDEL is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation; either
% version 2 of the License, or (at your option) any later version.
if ~exist('params', 'var')
	error('error: params not available');
end
if ~isa(params, 'struct')
	error('error: params is not a struct');
end
if ~isfield(params, 'Tcold1')
	error('error: symbolic parameter Tcold1 not defined in params structure');
end
if ~isfield(params, 'Tcold2')
	error('error: symbolic parameter Tcold2 not defined in params structure');
end
if ~isfield(params, 'Thot1')
	error('error: symbolic parameter Thot1 not defined in params structure');
end
if ~isfield(params, 'Thot2')
	error('error: symbolic parameter Thot2 not defined in params structure');
end
if ~isfield(params, 'Ts')
	error('error: symbolic parameter Ts not defined in params structure');
end
if ~isfield(params, 'Uc')
	error('error: symbolic parameter Uc not defined in params structure');
end
if ~isfield(params, 'Uh')
	error('error: symbolic parameter Uh not defined in params structure');
end
if ~isfield(params, 'alpha1')
	error('error: symbolic parameter alpha1 not defined in params structure');
end
if ~isfield(params, 'alpha2')
	error('error: symbolic parameter alpha2 not defined in params structure');
end
if ~isfield(params, 'k1')
	error('error: symbolic parameter k1 not defined in params structure');
end
if ~isfield(params, 'k2')
	error('error: symbolic parameter k2 not defined in params structure');
end
if ~exist('x', 'var')
	error('error:  current state x not supplied');
end
x=x(:);
if ~all (size(x)==[2 1])
	error('error: state vector has wrong dimension');
end
if ~exist('u', 'var')
	error('error: input u not supplied');
end
u=u(:);
if ~all (size(u)==[1 1])
	error('error: input vector has wrong dimension');
end

d = zeros(6, 1);
z = zeros(2, 1);
xn = zeros(2, 1);
y = zeros(2, 1);

if (x(1) < -10) | (x(1) > 50)
	error('variable T1 is out of bounds');
end
if (x(2) < -10) | (x(2) > 50)
	error('variable T2 is out of bounds');
end
if (u(1) < -10) | (u(1) > 50)
	error('variable Tamb is out of bounds');
end

% cold1 = T1 <= Tcold1;
within((x(1)) - (params.Tcold1), min((-10) + (-params.Tcold1), (50) + (-params.Tcold1)), max((50) + (-params.Tcold1), (-10) + (-params.Tcold1)), 28);
if (x(1)) - (params.Tcold1) <= 0
	d(3) = 1;
else
	d(3) = 0;
end

% cold2 = T2 <= Tcold2;
within((x(2)) - (params.Tcold2), min((-10) + (-params.Tcold2), (50) + (-params.Tcold2)), max((50) + (-params.Tcold2), (-10) + (-params.Tcold2)), 29);
if (x(2)) - (params.Tcold2) <= 0
	d(4) = 1;
else
	d(4) = 0;
end

% hot1 = T1 >= Thot1;
within((params.Thot1) - (x(1)), min(((0) + (params.Thot1)) + (-50), ((0) + (params.Thot1)) + (10)), max(((0) + (params.Thot1)) + (10), ((0) + (params.Thot1)) + (-50)), 26);
if (params.Thot1) - (x(1)) <= 0
	d(1) = 1;
else
	d(1) = 0;
end

% hot2 = T2 >= Thot2;
within((params.Thot2) - (x(2)), min(((0) + (params.Thot2)) + (-50), ((0) + (params.Thot2)) + (10)), max(((0) + (params.Thot2)) + (10), ((0) + (params.Thot2)) + (-50)), 27);
if (params.Thot2) - (x(2)) <= 0
	d(2) = 1;
else
	d(2) = 0;
end

% uhot = {IF cold1 | (cold2 & ~hot1) THEN Uh ELSE 0};
d(5) = (d(3)) | ((d(4)) & (~d(1)));

% ucold = {IF hot1 | (hot2 & ~cold1) THEN Uc ELSE 0};
d(6) = (d(1)) | ((d(2)) & (~d(3)));

% ucold = {IF hot1 | (hot2 & ~cold1) THEN Uc ELSE 0};
if d(6)
	within(params.Uc, min((0) + (params.Uc), (0) + (params.Uc)), max((0) + (params.Uc), (0) + (params.Uc)), 32);
	z(2) = params.Uc;
else
	within(0, 0, 0, 32);
	z(2) = 0;
end

% uhot = {IF cold1 | (cold2 & ~hot1) THEN Uh ELSE 0};
if d(5)
	within(params.Uh, min((0) + (params.Uh), (0) + (params.Uh)), max((0) + (params.Uh), (0) + (params.Uh)), 31);
	z(1) = params.Uh;
else
	within(0, 0, 0, 31);
	z(1) = 0;
end

% T1 = T1 + Ts * (-alpha1 * (T1 - Tamb) + k1 * (uhot - ucold));
xn(1) = (x(1)) + ((params.Ts) * (((-params.alpha1) * ((x(1)) - (u(1)))) + ((params.k1) * ((z(1)) - (z(2))))));

% T2 = T2 + Ts * (-alpha2 * (T2 - Tamb) + k2 * (uhot - ucold));
xn(2) = (x(2)) + ((params.Ts) * (((-params.alpha2) * ((x(2)) - (u(1)))) + ((params.k2) * ((z(1)) - (z(2))))));

% y1 = T1;
y(1) = x(1);

% y2 = T2;
y(2) = x(2);

xn=xn(:);
y=y(:);
z=z(:);
d=d(:);


function within(x, lo, hi, line)
 if x<lo | x>hi 
 error(['bounds violated at line ', num2str(line), ' in the hysdel source']); 
 end
