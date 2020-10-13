function pt = arc(ri, ro, a1, a2, N, equal_step)
% pt = arc(ri, ro, a1, a2, N, equal_step)
% generates a two-column coordinate vector pt, with the first column being
% the x-coordinates, second - y, of an arc with the inner radius, ri, and
% outer - ro. Additionally, ri and ro can be vectors of length 2 to prodice
% an elliptical boundary with the first element giving the semi axis along
% the x-direction and the second -- y-direction.
% Calls the ellipse function twice to generate each boundary and then
% stiches them together. 
% a1 and a2 - angles of the arc element (in radians)
% N - vector of length 2 with the numbers of points for the inner and outer 
% arcs (can be a scalar for an equal number). Default: 16 or
% abs(a1-a2)/pi*90 (every 2 degrees), whatever's larger.
% equal_step -- t/f interpolate the points on the curve to have an equal
% step accoring to the distance on the curved path and not angle. May be time 
% consuming. Default: false
%
% See also: ellipse, interparc

if nargin < 6, equal_step = false; end
if nargin < 5 || isempty(N)
    N = max(16, abs(a1-a2)/pi*90);
end
c = [0 0]; % coordinate center
pt = [ellipse(ri(1), ri(end), c, linspace(a1, a2, N(1)), 0, equal_step); ...
    flipud(ellipse(ro(1), ro(end), c, linspace(a1, a2, N(end)), 0, equal_step))];
% close the path
pt = [pt; pt(1,:)];
