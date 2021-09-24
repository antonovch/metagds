function L = duplicate_gds(L, N, R_trans)
% duplicates the TOP structure in the gds_library L onto a rectangular
% N(1) x N(end) grid with offset R_trans.

    % get all structs from the library
    st = get(L);
    % find index of the top layer
    istop = strcmp(cellfun(@sname,st,'UniformOutput',false),'TOP');
    % create cell array for all copies
    el_trans = cell(N);
    % get all elements of the top struct
    el_top = get(st{istop});
    for i = 1:size(el_trans,1)
        for j = 1:size(el_trans,2)
            el_trans{i,j} = el_top;
            for k = 1:length(el_top)
                el_trans{i,j}{k} = set(el_trans{i,j}{k}, ...
                    'xy', get(el_trans{i,j}{k}).xy + ([i,j]-1).*R_trans);
            end
        end
    end
    top = gds_structure('TOP', [el_trans{:}]); 
    L = gds_library('META.DB', 'uunit',get(L,'uunit'), ...
        'dbunit',get(L,'dbunit'), [st(~istop), {top}]);
end