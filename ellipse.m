function pt = ellipse(a, b, c, t, angle, equal_step)
% pt = ellipse(a, b, c, t, angle, equal_step)
% Generates a two-column coordinate vector pt, with the first column being
% the x-coordinates, second - y, of the ellipse with the radius along
% x-direction being given by a, and y-direction - by b. If a==b (the case when
% b is not given explicitely), the shape is a circle. 
% Additional optional arguments:
% c -- 2-element vector of x and y corrdinates of the center. Default: [0, 0]
% t -- vector (or a scalar giving its size) of angle values between 0 and
% 2*pi for the point definition. Default: 64.
% angle -- rotation angle in degrees. Default: 0.
% equal_step -- t/f interpolate the points on the curve to have an equal
% step accoring to the distance on the curved path and not angle (if t has
% an equal step). May be time consuming. Default: false
%
% See also: interparc
    if nargin < 6; equal_step = false; end
    if nargin < 5; angle = 0; end
    if nargin < 4
        t = linspace(0,2*pi,64).'; 
    elseif isscalar(t)
        t = linspace(0,2*pi,t).'; 
    else
        t = t(:);
    end
    if nargin < 3; c = [0 0]; else, c = c(:).'; end
    if nargin < 2; b = a; end
    if a == 0 || b == 0; pt = c; return; end
    pt = [a*cos(t), b*sin(t)];
    if angle ~= 0
        pt = [cosd(angle)*pt(:,1) - sind(angle)*pt(:,2),...
              sind(angle)*pt(:,1) + cosd(angle)*pt(:,2)];
    end
    pt = pt + c;
    if equal_step
        pt = interparc(length(t), pt(:,1), pt(:,2), 'pchip');
    end
end