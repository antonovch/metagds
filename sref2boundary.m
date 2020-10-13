function el = sref2boundary(sref, st)
% el = sref2boundary(sref, st)
% Transforms sref into boundary element. For this, we take the first
% argument, sref, which is an sref type gds_element instance and gives the 
% transofrmation (translation, rotation) that we need to perform on some 
% structure, which, in turn, has to passed as a second argument, st. 
% Thus, st is an instance of gds_structure. It also can be a cell array of 
% different structures and we'll choose the one referenced to autmatically.
% This strucure gives the original path to be transofrmed and recorded as a
% new bounary-type gds_element el, which is the output of this function.
% See also: gds_element, gds_structure
    sref = struct(sref);
    % is cell array given, find the one that sref references
    if iscell(st)
        snames = cellfun(@sname, st, 'UniformOutput', false);
        st = st{strcmp(snames, sref.data.sname)};
    end
    % extract gds_element(s) from the gds_structure
    el = struct(st).el; % cell array even if single element
    for ii = 1:length(el)
        el{ii} = bndtrans(el{ii}, sref);
    end
    % the translations in sref can be a two column vector, where the number
    % of row means a new element, so el will be a cell array of cell arrays
    % that has to be flatten out.
    el = [el{:}];
end
    
function cel = bndtrans(el, sref)
% subroutine that makes the actual translation. Two parameters are at play:
% strans.angle, which is a rotation of a each element, and data.xy, which
% is a two column vector, where each row is a (x,y) coordinate pair
% translation for an individual element, i.e. that's how you specify
% multiple new elements that are translations of a single fundumental one
    el = struct(el);
    % data.xy is a cell array that contains different paths that make up
    % this gds_element instance. Rotate each one by angle in radians
    rot = @(pt,ang)[cos(ang)*pt(:,1) - sin(ang)*pt(:,2),...
                        sin(ang)*pt(:,1) + cos(ang)*pt(:,2)];
    el.data.xy = cellfun(@(x)rot(x, sref.data.strans.angle), el.data.xy, 'UniformOutput', false);
    cel = cell(1,size(sref.data.xy,1));
    for ii = 1:size(sref.data.xy,1)
        cel{ii} = gds_element('boundary', 'xy', cellfun(@(x) x+sref.data.xy(ii,:), el.data.xy, 'UniformOutput', false),...
            'layer', el.data.layer);
    end
end
