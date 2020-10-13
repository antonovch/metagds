function [tls, opts] = gen_gds_structure(opts, I, R)
% [tls, opts] = gen_gds_structure(opts, I, R)
% Generates the basic structures to be combined to make the final gds file.
% Inputs:
%   opts - options struct from genGDS function
%   I - vector with numerical indicies of the stuctures to generate. By
%   default, all the geometries in opts.geom
%   R - current radius of cirvature
% Outputs:
%   tls - cell array of gds_structures (or empty vectors) for the
%   corresponding geometries and specifications from opts. It can be empty
%   when the structure has to be re-generated for each R
%   opts - same struct. Note that it might change, namely, in one of the
%   operation regimes, opts.geom elements can be substituted with a
%   gds_structure to be further adapted for a specific R.
%
% See also: genGDS, ellipse, arc, rectbox, gds_structure, gds_element, 
% eldiff, sref2boundary

    if nargin < 3; R = 0; end
    if nargin < 2 || isempty(I); I = 1:numel(opts.geom); end
    tls = cell(1,numel(I));
    for iii = 1:numel(I)
        ii = I(iii);
        % width of the ring = y coordinate of a unit cell (imagine a unit cell 
        % sitting at 12 o'clock, we're building rings from 3 o'clock, so PR
        % becomes length along x)
        PR = opts.period(ii,2); 
        PA = opts.period(ii,1); 
        Rmid = R + PR/2; % radius of the middle of the unit cell
        numSec = round(2*pi*Rmid/PA); % number of section in the ring (determines the angular spread of the arc elements)
        alpha = 2*pi/numSec; % angular width of the arc elemens
        % function handle for rotating a boundary
        rot = @(pt,ang)[cos(ang)*pt(:,1) - sin(ang)*pt(:,2),...
                        sin(ang)*pt(:,1) + cos(ang)*pt(:,2)];
        
        geom = opts.geom{ii};
        if iscell(geom)
            cls = class(geom{1});
        else
            cls = class(geom);
        end
        switch cls
            case 'struct'
                % filling missing field, most of interst are c and angle
                geom = arrayfun(@(s)updatestruct(struct('r',0,'c',[0 0],'fun','ellipse','angle',0), s),geom);
                xy = cell(1,numel(geom)); 
                if strcmp(opts.type, 'circ+')
                % creating cell arrays of gds_structures whose position will 
                % later be adjusted for a partucliar R (curvature) in the
                % 'cell' case of this switch.
                    % step 1: get the paths
                    for kk = 1:length(xy)
                        if ischar(geom(kk).fun)
                            switch geom(kk).fun
                                case 'ellipse'
                                    xy{kk} = ellipse(geom(kk).r(1), geom(kk).r(end), ...
                                        [0 0], opts.pts, geom(kk).angle, opts.uniform_boundary);
                                case 'rectbox'
                                    xy{kk} = rectbox(geom(kk).r(1), geom(kk).r(end), [0 0], geom(kk).angle);
                            end
                        else
                            xy{kk} = feval(geom(kk).fun, geom(kk).r{:}, [0 0]);
                        end
                    end
                    % step 2: make required parent/child substractions
                    xy = substract_boundaries(xy, opts.A{ii}, opts.n(ii));
                    % step 3: transform paths to gds_struct
                    for kk = 1:length(xy)
                        xy{kk} = gds_structure(['BLOCK_',num2str(ii),'_',num2str(kk)], ...
                                gds_element('boundary', 'xy', xy{kk}, 'layer', opts.layer));
                    end
                    opts.geom{ii} = [xy, geom];
                    continue % carry on, leaving tle/tls empty
                else
                    % creating the whole unit cell (inclusions are on a rectangular grid)
                    for kk = 1:length(xy)
                        if ischar(geom(kk).fun)
                            switch geom(kk).fun
                                case 'ellipse'
                                    xy{kk} = ellipse(geom(kk).r(1), geom(kk).r(end), ...
                                        [Rmid 0]+geom(kk).c, opts.pts, geom(kk).angle, opts.uniform_boundary);
                                case 'rectbox'
                                    xy{kk} = rectbox(geom(kk).r(1), geom(kk).r(end), ...
                                        [Rmid 0]+geom(kk).c, geom(kk).angle);
                            end
                        else
                            xy{kk} = feval(geom(kk).fun, geom(kk).r{:}, [Rmid 0]+geom(kk).c);
                        end
                    end
                    % making required parent/child substractions
                    xy = substract_boundaries(xy, opts.A{ii}, opts.n(ii));
                    tle = gds_element('boundary', 'xy', xy, 'layer', opts.layer);
                end
            case 'gds_structure'
                % uses previously pre-generated gds_sctructres for the
                % 'circ+' type and makes a unit cell out of them.
                % geom = cell array of gsd_structures, except for the last,
                % which is a stuct.
                assert(nargin > 2, 'Requires the radius as an argument to update the unit cell');
                s = geom{end};
                tle = cell(1,numel(geom)-1);
                for kk = 1:numel(tle)
                    % the additional circumference
                    circ_diff = 2*pi*s(kk).c(2); % 2*pi*(Rmid+center_y) - 2*pi*Rmid
                    % per section and halved because we move from the
                    % center into one direction
                    circ_diff = circ_diff/numSec/2;
                    % add scaling factor that depends on the distance from
                    % the central axis, plus it also makes sure that the angle
                    % has the correct sign, which is sign(center_x)*sign(center_y)
                    circ_diff = circ_diff*(s(kk).c(1)/PA);
                    % get the angle (ang*R=arc_len)
                    ang = circ_diff/(Rmid+s(kk).c(2));
                    tle{kk} = gds_element('sref', 'sname', sname(geom{kk}), ...
                        'xy',rot([Rmid 0]+s(kk).c, ang), 'strans',struct('angle',0));
                end
            otherwise % all numerical options
                % geom is either an image, a path, or cell array of path
                if ~iscell(geom) 
                    if size(geom, 2) ~= 2
                        % image given, need to extract boundary
                        xgrid = linspace(0, PA, size(geom,1)).' - PA/2; 
                        ygrid = linspace(0, PR, size(geom,2)).' - PR/2; 
                        % compute boundaries and adjacency matrix
                        [b,~,n,A] = bwboundaries(round(geom),8,'holes');
                        opts.n(ii) = n;
                        opts.A{ii} = A;
                        % transform from indicies to coordinates
                        for ib = 1:length(b)
                            b{ib} = [xgrid(b{ib}(:,1)), ygrid(b{ib}(:,2))];
                        end 
                        geom = b(:).';
                    else
                        % single path, make it a cell for compatibility
                        geom = {geom}; 
                    end
                end
                % substract children form parents 
                geom = substract_boundaries(geom, opts.A{ii}, opts.n(ii));
                
                if strcmp(opts.type, 'circ+') 
                    s = struct('c',cellfun(@(x)mean(x,1), geom,'UniformOutput',0));
                    for kk = 1:length(geom)
                        geom{kk} = gds_structure(['BLOCK_',num2str(ii),'_',num2str(kk)], ...
                            gds_element('boundary', 'xy', geom{kk}, 'layer', opts.layer));
                    end
                    opts.geom{ii} = [geom, s];
                    continue % carry on, leaving tle/tls empty
                else
                    % Rmid here should be normally equal to PR/2
                    if contains(opts.type, 'cyl')
                        geom = cellfun(@(x) x+[PA/2 0], geom, 'Uniformoutput', false); 
                    else
                        geom = cellfun(@(x) x+[Rmid 0], geom, 'Uniformoutput', false);
                    end
                    tle = gds_element('boundary', 'xy', geom, 'layer', opts.layer);
                end
        end
        if (opts.ispost(ii) && strcmp(opts.shading, 'holes')) || (~opts.ispost(ii) && strcmp(opts.shading, 'posts'))
        % invert the shading of the unit cell element by creating a uniform
        % backround and then substructing tle from it. The complexity lies
        % in the fact that GDS goes not permit polygons with holes, i.e. to 
        % create a torus one need to connect outer with inner boundaries, 
        % like drawing without removing your hand. We do this with eldiff 
        % function. However, it's rather limited and might not work for all
        % cases. Things to try are to split the backgound into smaller tiles,
        % or use logical operation functions from the GDSII toolbox,
        % which might work funny. For one, this tiling is practically required.
        % And two, it sometimes throws error due to a hole, which is not
        % really there, so you could just suppress this error in the code
        % manually. If eldiff fails, the code with try GDSII function.
        % There is no option at the moment to start with the latter in the
        % first place (in case eldiff does not fail, but produces incorrect
        % result).
        
            % if it's a circlular pattern, we have to make arcs, instead of
            % rectangles, and their size depends on the radius
            N = opts.bgd_tiling;
            if isscalar(N); N = N.*[1 1]; end
            if strcmp(opts.type, 'circ') || strcmp(opts.type, 'circ+')
                % in circ mode, an appropriate full arc element is to be
                % generated only when R is given in the process of executing 
                % genGDS's while loop
                if nargin < 3
                    continue
                end
                W = PR/N(1); % width of each sub section
                % inner radii
                R_inner = Rmid + (-N(1)/2:N(1)/2-1) * W;
                angle_step = alpha/N(2);
                angle_start = angle_step * (-N(2)/2:N(2)/2-1);
                bg = cell(length(R_inner),length(angle_start));
                for jj = 1:length(R_inner)
                    for kk = 1:length(angle_start)
                        bg{jj,kk} = arc(R_inner(jj), R_inner(jj)+W,...
                            angle_start(kk), angle_start(kk)+angle_step,...
                            ceil(opts.pts/2), opts.uniform_boundary);
                    end
                end
                bg = gds_element('boundary','xy',bg(:).','layer', opts.layer);
                
                %%% old code using gdsii_arc function 
%                 bg = gds_element('boundary','xy',{},'layer', opts.layer);
%                 for jj = -N(2)/2+1:N(2)/2
%                     for kk = -N(1)/2:N(1)/2-1
%                         arc_struct = struct('r', Rmid+kk*PR/N(1),'c',[0 0],...
%                         'a1',alpha/N(2)*(jj-1),'a2',alpha/N(2)*jj,...
%                         'w',PR/N(1),'e',1e-6,'pap',0);
%                         bg = bg + gdsii_arc(arc_struct, opts.layer);
%                     end
%                 end
            else
                [X,Y] = ndgrid(0:N(1)-1,-N(2)/2:N(2)/2-1);
                % size = 5x2xprod(N)
                xy = ([0 0; 1 0; 1 1; 0 1; 0 0] + reshape([X(:), Y(:)].', 1, 2, [])).*[PA PR]./N;
                % break into cell array along the 3rd dim
                xy = squeeze(mat2cell(xy,size(xy,1),2,ones(1,size(xy,3))));
                bg = gds_element('boundary', 'xy',xy, 'layer', opts.layer);
            end
            if iscell(tle) % array of sref elements
                for ll = 1:length(tle)
                    el = sref2boundary(tle{ll}, opts.geom{ii});
                    for tt = 1:length(el)
                        try
                            bg = eldiff(bg, el{tt});
                        catch
                            bg = poly_bool(bg, el{tt}, 'notb', 'layer', opts.layer);
                        end
                    end
                end
                tle = bg;
            else % single boundary element
                try
                    tle = eldiff(bg, tle);
                catch
                    tle = poly_bool(bg, tle, 'notb', 'layer', opts.layer);
                end
            end
        end
        tls{iii} = gds_structure(['EL_', num2str(ii)], tle);
    end
end

function xy = substract_boundaries(xy, A, n)
    if nnz(A) == 0 || n == 0, return, end
    for ia = 1:n
        for ib = find(A(:,ia)).'
            xy{ia} = eldiff(xy{ia},xy{ib});
            xy(ib) = {[]}; % distinct form xy(ib) = [], 
            % works like xy{ib} = [], but also tolerates empty ib
        end
    end
    xy(cellfun(@isempty, xy)) = [];
end
