function pt = rectbox(a, b, c, angle)
% pt = rectbox(a, b, c, angle)
% All arguments, except for "a" are optional.
% generates a two-column coordinate vector pt, with the first column being
% the x-coordinates, second - y, of the square/rectangle box with the width
% (x-direction) given by a, and height (y-direction) - by b. 
% size(pt,1)==5 if neither a or b are zero, otherwise - 3.
% c -- 2-element vector of x and y corrdinates of the center. Default: [0, 0]
% angle -- rotation angle in degrees. Default: 0.
%
% Equivalent to calling ellipse(a,b,c,(1:2:9)*pi/4,angle)*sqrt(2) up to a
% rounding error.
%
% See also: ellipse

    if nargin < 4; angle = 0; end
    if nargin < 3; c = [0 0]; end
    if nargin < 2; b = a; end
    if a == 0 || b == 0
        pt = c + [a, b]./[-2; 2];
    else
        pt = [a, b].*[1 1; -1 1; -1 -1; 1 -1]/2 + c;
    end
    
    if angle ~= 0
        pt = [cosd(angle)*pt(:,1) - sind(angle)*pt(:,2),...
              sind(angle)*pt(:,1) + cosd(angle)*pt(:,2)];
    end
    % close the path
    pt = [pt; pt(1,:)];
end
