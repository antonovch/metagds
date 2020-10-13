function dif = eldiff(el1, xy2)
% An implementation of GDSII toolbox's poly_bool function with 'notb'
% option aka difference between the two elements. It's pretty simple, so
% there are probably a lot of cases where it could fail, but if you use
% ellipse function to generate your contors, it seems OK. Still have to do
% with the gds_element class non-sence. Tested with the older version of
% the toolbox. 
% el1 - a scalar instance of gds_element or a path
% xy2 - cell array of paths (two column vector of x and y coordinates), 
% gds_elements, or a single gds_element instance.
% dif - the result of taking the difference. Same class as el1.
% Example simple use (poly_bool fails because 'a hole was created'):
% bg = gds_element('boundary', 'xy', ellipse(20), 'layer', 6);
% el1 = gds_element('boundary', 'xy', ellipse(5,6,[10,6]), 'layer', 6);
% el2 = gds_element('boundary', 'xy', ellipse(2,2,[-10,6]), 'layer', 6);
% bg=eldiff(bg,el1);
% tst=eldiff(bg,el2);
% write_gds_library(gds_library('TEST.DB','uunit',1e-6, 'dbunit',1e-11, ...
%     gds_structure('TEST', {tst}), 'layer',6), '!test.gds');

    if nargin < 2 || isempty(xy2), dif = el1; return; end
    gdsel = false;
    if isa(el1, 'gds_element')
        xy1 = struct(el1).data.xy;
        assert(length(xy1) == 1, 'Only for single path gds_element')
        xy1 = xy1{1};
        gdsel = true;
    elseif iscell(el1)
        assert(length(el1) == 1, 'Only for single path gds_element')
        xy1 = el1{1};
    else
        xy1 = el1;
    end
    
    switch class(xy2)
        case 'gds_element'
            xy2 = struct(xy2).data.xy;
        case 'cell'
            if isa(xy2{1}, 'gds_element')
                xy2 = cellfun(@(x) struct(x).data.xy, xy2, 'UniformOutput', false);
                xy2 = [xy2{:}];
            end
        case 'double'
            xy2 = {xy2};
        otherwise
            error('Boundary of class %s is not supported', class(xy2))
    end
    
    % make sure they're both counter-clockwise. They can be whatever, we
    % just need to know to make the correct path, and ellipse() outputs
    % them anti-clockwise.
    if ispolycw(xy1(:,1), xy1(:,2)), xy1 = flipud(xy1); end
    iscw = cellfun(@(xy) ispolycw(xy(:,1), xy(:,2)), xy2);
    xy2(iscw) = cellfun(@flipud, xy2(iscw), 'UniformOutput', false);
    
    % find two closest points between xy1 and each of xy2 paths
    i1 = zeros(1,length(xy2));
    for ii = 1:length(xy2)
        dist = sqrt((xy1(:,1)-xy2{ii}(:,1).').^2 + (xy1(:,2)-xy2{ii}(:,2).').^2);
        [~, I] = min(dist(:));
        [i1(ii),i2] = ind2sub([size(xy1,1) size(xy2{ii},2)], I);
        xy2{ii} = circshift(xy2{ii},-i2+1);
    end
    [i1, is] = sort(i1);
    xy2 = xy2(is);
    i1 = [i1 size(xy1,1)];
    % dif contains parts of the path that is the difference between the elements, 
    % later to be contatenated
    dif = cell(1,length(i1)); 
    dif{1} = xy1(1:i1(1),:);
    for ii = 1:length(i1)-1
        xy = [xy2{ii}; xy2{ii}(1,:)];
        dif{ii+1} = [flipud(xy); xy1(i1(ii):i1(ii+1),:)];
    end
    dif = cat(1, dif{:});
    if gdsel % if gds_element was given as xy1, return a gds_element too
        dif = gds_element('boundary', 'xy', dif, 'layer', el1.layer);
    end
end