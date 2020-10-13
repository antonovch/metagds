function L = genGDS(varargin)
% L = genGDS(phase_data, varargin)
% L = genGDS(varargin)
% Generates a GDS file with metasurface/metagrating geometry.
% Inputs:
% - phase_data -- phase distribution map as a function of a
%   single coordinate (axial/cylindrical symmetry). 
%   Can be one of:
%   1/ struct with two fields: R and phase, both vectors of 
%      equal size containing the corresponding data;
%   2/ two column matrix with R and phase values; 
%   3/ function_handle that takes R and returns phase. At
%   the moment, phase map can vary only radially (circular, 
%   cylindrical). 
% - varargin -- Name/Value pairs with the following options 
%               (deafult values in [brackets]):
%   
%   * type, char array, ['circ2cart']
%       The type of the surface pattern. Available options are 'circ', 
%       circular metasurface on a polar grid, using the METAC algorithm;
%       'circ+', same as 'circ', but tries to adjust the positions of  
%       elements of a given unit cell structure to account for the distortions  
%       due to the curvature, might cause errors; 'circ2cart', circular 
%       pattern on a Cartesian grid; 'cyl' and 'cyl2cart' cylindrical 
%       metasurfaces with different internal structure, 'cyl' is generally
%       preferred.
%   * shading, 'holes' or 'posts', ['posts']
%       specifies whether the shaded region of the GDS file should 
%       correspond to the holes, or posts (areas occupied by the material) 
%       in the design.
%   * fname, char array, optional
%       name of the output file. Adding the format (`.gds') is optional. 
%       Prepend the name with `!' to overwrite a potentially existing file 
%       of the same name. If empty, no file is created.
% 
%   * geom, cell array, required
%       a cell array (or numerical matrix/struct for a single unit cell), 
%       where each cell defines a unique meta-atom geometry. There are three
%       ways to defile a geometry: (i) numerical or logical array representing 
%       a binary image, (ii) two-column matrix of x and y coordinates tracing 
%       a polygon, (iii) struct with the following fields:
%           - fun, name of the function to use to get the boundary, either 
%             'ellipse' or 'rectbox' to use the corresponding functions, 
%              or a custom function handle or name. 
%           - r, two element vector with semi-axes radii if fun is 'ellipse',
%              or width and height pair if it's 'rectbox', or a cell array
%              of the first N-1 arguments of the custom functions, out of N.
%           - c, coordinates of the element's center within the unit cell.
%              If a custom function is given, has to be the final argument.
%           - angle, angle of rotation (in degrees) around c.
%     
%    * period, numerical array, depends
%       generally, a two column vector of periods in x and y directions, 
%       with rows corresponding to respective elements in geom, and required, 
%       but for image inputs can be omitted -- their sizes will be used instead, 
%       i.e. one pixel is one uunit large. Can be a single row if all basic 
%       unit cells have the same periods, and single column for a square unit cell.
%     
%   * phi0, numerical array, depends
%       numerical array of characteristic phase responses in radians of unit 
%       cells in geom. Optional only in test mode.
%     
%   * uniform_boundary, logical, [false]
%       defines the equal_step parameter of ellipse function for the case 
%       when geom is defined as a struct and geom.fun is 'ellipse'). 
%       It uses the interparc} function in order to select points on a 
%       curved surface with equal separation along the trajectory. 
%       The function can be downloaded from https://www.mathworks.com/matlabcentral/fileexchange/34874-interparc
%     
%   * A, cell array, optional
%       cell array of the same size as geom (or just a matrix for a single 
%       meta-atom) of either full or sparse square (usually logical) matrices 
%       that represent the parent-child dependencies between object boundaries 
%       and hole boundaries. For a unit cell defined as either a cell array 
%       of paths, or an array of struct, A gives an adjacency matrix between 
%       the boundaries. For unit cells that are defined through an image, 
%       A is calculated automatically with bwboundaries. Therefore, its 
%       definition follows the one in the output of bwboundaries. If the 
%       matrix has a non-zero element (true) at the position (i,j), it 
%       means that ith boundary (child) will be subtracted from jth (parent). 
%       If there is no parent/child elements for a given unit cell, its A 
%       can be either empty or all zeros.
%     
%   * n, numerical/cell array, optional
%       array (size equal to geom's) of integer numbers of objects in the 
%       unit cells. This is bwboundaries's definition that does not make much 
%       sense here. What to do with it, depends on the way geom is defined:
%       1/ geom is a binary image: n is not required, as it's calculated 
%       automatically using bwboundaries, along with A;
%       2/ It's a struct and A is constructed manually: in this case, 
%       n would be equal to the number of columns (or rows, since A is a 
%       square matrix) in the corresponding A, as the user specifies all 
%       the subtractions needed. In this case, the user does not have to 
%       explicitly provide n.
%       3/ It's a cell array of boundaries: if A was defined manually, 
%       n similarly can be omitted, if, however, it was obtained using 
%       bwboundaries, n has to be fetched from the eponymous output variable, 
%       since, in general, size(A,2) would be larger than n, as 
%       bwboundaries first puts objects (parents), then holes (children), 
%       and we only need to loop over the objects.
%   * uunit, numerical scalar, [1e-6]
%       user units as a fraction of a meter. Default -- microns.
%   * dbunit, numerical scalar, [1e-11]
%       database units as a fraction of a meter, used for storing the data. 
%       In other words, it determines the resolution. Keep in mind, that 
%       GDS stores the positions as a 4-byte integer, so there's a trade-off 
%       between the resolution and the maximum size (largest coordinate value) 
%       of the structure. The number of significant digits after the decimal 
%       marker in a length specified in user units is given by log10(uunit/dbunit).
%   * Rmin, numerical scalar, [0]
%       inner radius of the metasurface. For a circular metasurface it means 
%       that we have a void in the center of radius Rmin, which is 
%       automatically shaded if shading parameter is 'holes'. This might be 
%       useful to avoid regions of extreme curvature. For a cylindrical 
%       metasurface, it shows a shift of coordinates in x-direction, which 
%       is not useful and should be left 0.
%   * Rmax, numerical scalar or two element vector, depends
%     	outer radius of the metasurface. For a cylindrical (rectangular) 
%       metasurface, it is a two-element vector of the upper-right corner 
%       of the device. Optional only in test mode.
%   * ispost}, logical, [true]
%       specifies whether the provided boundary in geom encircles areas 
%       of the material, or, if an image is given, whether non-zero pixels 
%       mark the (higher-index) material.
%   * pts, numerical scalar, [64]
%   	number of points per single boundary generated by the ellipse 
%       function (if it's used at all), i.e. its t argument.
%   * bgd_tiling, numerical array, [1]
%       size of grid to break the 'background' of a unit cell into. 
%       Can be a scalar for a square grid. If we need to produce a "full 
%       unit cell", gen_gds_structure} takes the 'background', which 
%       corresponds to the metasurface unit cell, which is either a rectangle 
%       or an arc, with filling fraction equal to one, and subtracts from 
%       it all subelements whose boundaries are provided.
%   * layer, numerical scalar, [1]
%       layer number in the GDS file. Not really important.
%   * label, struct, optional
%       parameters of the label to be created. A struct with the following fields:
%         - string, text to be written;
%         - height, maximum letter height in user units (uunit);
%         - max_feature, largest feature i.e. the square diameter for the 
%           checkerboard pattern;
%         - location, position on the label, either a two-coordinate vector 
%           of the lower left corner, or one of the following options:
%               • 'ne': north-east, upper right corner,
%               • 'nw': north-west, upper left corner,
%               • 'out': the text ends in the upper right corner right 
%                 outside the metasurface region for cylindrical metasurfaces 
%                 (rectangular device) and inside for circular (corners are empty),
%               • 'center': text hovers above the middle point of the device.
%           Optional, 'ne'} by default;
%         - mirror, boolean, whether to position a copy the text upside 
%           down on the opposite side of the device. Optional, true by default;
%         - angle, angle to rotate the text in degrees. Optional, 0 by default.
%
% Optputs:
% - L -- gds_library object, which is written on disk if fname is provided.
%
% See metagds_doc.pdf for more info and gds_examples.mlx 
% for some examples of using this function.
% 
% Below are the functions that are required if not every time, 
% then at least for some feature of this function.
%
% See also: gen_gds_structure, sref2boundary, eldiff, updatestruct, ellipse, 
% arc, interparc

    wb = waitbar(0, 'Starting genGDS...');
    opts = struct('geom', [], 'phi0', [], 'period', [], 'ispost', true, ...
        'Rmin', 0, 'Rmax', [], 'type', 'circ2cart', 'fname', '', ...
        'A', [], 'n', [], 'pts', 64, 'label',...
        struct('string','', 'height', [], 'max_feature', [], 'angle', 0, ...
        'mirror', true, 'location', 'ne'),...
        'shading', 'posts', 'uunit', 1e-6, 'dbunit', 1e-11,...
        'layer', 1, 'bgd_tiling', 1, 'uniform_boundary', false); % default options
    if ~ischar(varargin{1}) && ~isstruct(varargin{1})
        phase_data = varargin{1};
        varargin(1) = [];
    else
        phase_data = [];
    end
    if isstruct(varargin{1}) % struct of options is given
        varargin = varargin{1};
    else
        assert(mod(numel(varargin),2)==0, 'Uneven number of input arguments')
        % this trick compensates the cell unpacking that struct performs, e.g.
        % struct('a',{1}) will register number 1 in the field 'a', not cell
        % arary {1}
        I = cellfun(@iscell, varargin);
        varargin(I) = cellfun(@(x){x}, varargin(I), 'UniformOutput', false);
    end
    opts = updatestruct(opts, varargin); % get user-supplied options
    gdsii_units(opts.uunit, opts.dbunit) % set global vars, not cool, but that's how the toolbox is
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% some checks and defaults 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~iscell(opts.geom), opts.geom = {opts.geom}; end
    M = numel(opts.geom);
    if isempty(opts.phi0) && isempty(phase_data)
        % testing mode
        opts.phi0 = linspace(0,pi,M);
    elseif numel(opts.phi0) ~= M
        error('Length of phi0 vector is not equal to the one of geom')
    end
    [opts.phi0, I] = sort(mod(opts.phi0, 2*pi));
    opts.geom = opts.geom(I);
    if numel(opts.ispost) == 1
        opts.ispost = repmat(opts.ispost, 1, M);
    else
        opts.ispost = opts.ispost(I);
    end
    if isempty(opts.period)
        opts.period = cellfun(@size, opts.geom, 'UniformOutput', false);
        opts.period = cat(1, opts.period{:});
        if any(opts.period(:) <= 2)
            % arrays of structs/boundaries have a dimension of 1, single
            % boundary - of 2
            error('For non-image inputs, period is a required argument.')
        end
    elseif size(opts.period,2) == 1
        opts.period = repmat(opts.period, 1, 2);
    end
    if size(opts.period,1) == 1
        opts.period = repmat(opts.period, M, 1);
    else
        opts.period = opts.period(I,:);
    end
    if ~iscell(opts.A), opts.A = {opts.A}; end
    if numel(opts.A) < M
        if ~all(cellfun(@isempty,opts.A))
            warning(['Number of adjacency matrices is smaller than that of ';
                'the unit cells. Will be applied to first %g cells.'], numel(opts.A))
        end
        opts.A = [opts.A, cell(1,M-numel(opts.A))];
    elseif numel(opts.A) > M
        error('Too many cells in "A".')
    end

    if isempty(opts.n)
        opts.n = cellfun(@length, opts.A);
    else
        if iscell(opts.n), opts.n = [opts.n{:}]; end
        if numel(opts.n) < M
            warning(['Number of "n" parameters is smaller than that of ';
            'the unit cells. Will be applied to first %g cells.'], numel(opts.n))
            opts.n = [opts.n, zeros(1,M-opts.n)];
        elseif numel(opts.n) > M
            error('Too many entries for "n".')
        end
        assert(all(opts.n <= cellfun(@numel, opts.geom)), ...
            'Number of objects n is larger than the actual number in one of the geom cells.')
    end
    if isempty(phase_data)
        if isempty(opts.Rmax)
            opts.Rmax = opts.Rmin + sum(opts.period(:,2));
            if contains(opts.type,'cyl')
                opts.Rmax(2) = max(opts.period(:,1));
            end
        end
        if M == 1 % single geom
            phase_data = @(R) repmat(opts.phi0, size(R));
        else
            phase_data = @(R) interp1(linspace(opts.Rmin(1), opts.Rmax(1), M), ...
                opts.phi0, R, 'spline');
        end
    else
        switch class(phase_data)
            case 'double'
                if isempty(opts.Rmax); opts.Rmax = phase_data(end, 1); end
                phase_data = @(R) interp1(phase_data(:,1), phase_data(:,2), R, 'spline');
            case 'struct'
                if isempty(opts.Rmax); opts.Rmax = phase_data.R(end); end
                phase_data = @(R) interp1(phase_data.R, phase_data.phase, R, 'spline');
            case 'function_handle'
                assert(~isempty(opts.Rmax), 'Rmax has to be provided if phase_data is a function_handle')
            otherwise
                error('Wrong phase_data variable')
        end
    end
    if contains(opts.type,'cyl') && isscalar(opts.Rmax) 
        opts.Rmax = [1 1]*opts.Rmax;
    end
    % function that rotates coordinates in pt by the angle ang in degrees
    rotd = @(pt,ang)[cosd(ang)*pt(:,1) - sind(ang)*pt(:,2),...
                     sind(ang)*pt(:,1) + cosd(ang)*pt(:,2)];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Print approximate number of unit cells in the device 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mpr = mean(opts.period,1);
    switch opts.type
        case {'circ', 'circ+'}
            Ntot = sum(round(2*pi*(opts.Rmin:mpr(2):opts.Rmax)/mpr(1)));
        case 'circ2cart'
            Ntot = round(pi*(opts.Rmax.^2-opts.Rmin.^2)/prod(mpr));
        otherwise % cyl and cyl2cart
            Ntot = round(prod(opts.Rmax)/prod(mpr));
    end
    nums = [mod(fix(Ntot/1e6), 1e3),mod(fix(Ntot/1e3), 1e3),mod(Ntot, 1e3)];
    nums(nums == 0) = [];
    fprintf(['Total number of individual tiles (approximately) = ',...
        repmat('%.0f ',1,length(nums)),'\n'], nums);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Pre-calculation of base elements
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tic
    maxNumRings = ceil((opts.Rmax(1)-opts.Rmin)/min(opts.period(:,2)));
    rings=cell(1,maxNumRings+1); % +1 to accomodate the center element
    % rings is a cell array of references to sections from allsects, which
    % are like jigsaw pieces for the whole structure, each with an sref (in
    % rings), potantially with a translation, that we bring together at the
    % "TOP" gds cell.
    if opts.Rmin ~= 0 && contains(opts.type, 'circ') && strcmp(opts.shading, 'holes') % incert a hole in the middle
        rings{1} = {gds_element('boundary',...
            'xy',ellipse(opts.Rmin,opts.Rmin,[0 0],opts.pts,0,opts.uniform_boundary), 'layer', opts.layer)};
    end
    ir = 1; % iteration index for the while loop that allows to incorporate 
    % different periods (widths of rings) for 'circ' and 'cyl' options
    PR = min(opts.period(:,2)); % current period in the radial direction (y for Cartesian)
    areaFraction_old = 0; % for the waitbar
    waitbar(0, wb, 'Pre-defining unit cells...');
    [tls0, opts] = gen_gds_structure(opts); % creates structures located in
    % the origin, to be either translated later, or re-created with the correct
    % curvature ('circ+' type). For this latter case, gen_gds_structure 
    % will modify opts.geom putting there a gds_structure, while keeping
    % tls0 empty. Then, gen_gds_structure will be called again with a particular 
    % R (curvature) and the position of unit cell's blocks will be adjusted. 
    % We still need to save those fundamental elements in the final file though.
    % They are cecorded in the allels cell array. We have to use cell arrays  
    % instead of, say, gds_struct arrays, because the toolbox does not do 
    % ()-referencing correctly due to poorly overloaded subsref method.  
    % Moreover, other methods expect either cell arrys or single instances.
    switch opts.type
        case 'circ'
            % sections in allsects for the top cell will be rings of
            % different radii and potentially width
            allels = [cell(1, maxNumRings+1), ... unit cells for each R
                      tls0]; %  basic unit cells (meta-atoms)
            allsects = cell(1,maxNumRings);
        case 'circ+'
            % sections in allsects for the top cell will be rings of
            % different radii and potentially width
            allels = [opts.geom{:}];
            I = cellfun(@(x) isa(x, 'gds_structure'), allels);
            allels = [cell(1, maxNumRings+1), ... unit cells for each R
                      allels(I)]; % consituent elements 
            allsects = cell(1,maxNumRings);
        case 'cyl'
            % sections in allsects for the top cell will be 1D stripes
            % along y of different unit cells of potentially different
            % width (according to the corresponding unit cell period)
            allsects = cell(size(tls0));
            for isec = 1:numel(allsects)
                PY = min(opts.period(isec,2));
                Y = (PY/2:PY:opts.Rmax(2)).';
                X = zeros(size(Y));
                % get sref's of a unit cell that form a stripe
                cur_ref = gds_element('sref', 'sname', sname(tls0{isec}), ...
                                'xy',[X Y], 'strans',struct('angle',0));
                % transform sref to structure
                allsects{isec} = gds_structure(['STRIPE_', num2str(isec)], cur_ref);
            end
            allels = tls0; % save basic unit cells
        case {'circ2cart', 'cyl2cart'}
            % sections in allsects for the top cell will be parts (not
            % necessarily continuous) of the whole structure composed of
            % the corresponding unit cell. Note that all periods have be
            % equal for a uniform Cartesian grid
            allels = tls0;
            [allsects, rings] = deal(cell(1,M));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Main loop that builds the device layout (top level cell)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    waitbar(0, wb, '0% done');
    switch opts.type
        case {'circ2cart', 'cyl2cart'}
            % Only for the cases where all unit cells periods are equal
            waitbar(0, wb, 'Getting the coordinates...');
            if strcmp(opts.type, 'circ2cart')
                [X, Y] = ndgrid(-(opts.Rmax-PR/2):PR:opts.Rmax-PR/2);
                R = sqrt(X.^2+Y.^2);
                I = R < opts.Rmax(1) & R > opts.Rmin;
                X = X(I); Y = Y(I); R = R(I);
            else % cyl2cart
                PA = min(opts.period(:,1));
                [X, Y] = ndgrid(opts.Rmin+PA/2:PA:(opts.Rmax(1)-PA/2),...
                                PR/2:PR:opts.Rmax(2));
                R = X;
            end
            [R, ~, ic] = unique(R);
            phases = mod(phase_data(R), 2*pi);
            [~, lvls] = min(abs(opts.phi0(:).' - phases(:)),[],2);
            lvls = lvls(ic);
            lvls = reshape(lvls, size(X));
            for ii = 1:M
                I = lvls == ii;
                if all(~I), continue, end
                copies = gds_element('sref', 'sname', sname(tls0{ii}), 'xy',[X(I)-PR/2 Y(I)], 'strans',struct('angle',0));
                allsects{ii} = gds_structure(['KIND_', num2str(ii)], copies);
                rings{ii} = gds_element('sref', 'sname', sname(allsects{ii}), 'xy',[0 0], 'strans', struct('angle',0)); 
                waitbar(ii/M, wb, sprintf('%.0f%% done', ii/M*100));
            end
        case {'circ', 'circ+'} % METAC for circular patterning
            R = opts.Rmin; % R is the inner radius of the current ring
            while R < opts.Rmax
                % generating values for phase each time, in case the ring sizes are unequal
                cur_lvl = get_el_ind(opts, phase_data(R+PR/2)); % see below for definition
                PR = opts.period(cur_lvl, 2); % update value of the width of the current ring
                Rmid = R + PR/2; % distance from the origin to the middle of the unit cell
                numSec = round(2*pi*Rmid/opts.period(cur_lvl, 1));
        %         PA = 2*pi*Rmid/numSec; % perdiod in angular (length of arc that goes through the middle of a unit cell)
            %     alpha = acos(1-P^2/2/R^2);
                alpha = rad2deg(2*pi/numSec); % angular span of the unit cell, as
                % it is no longer rectangular, but arc-shaped
                % flags is a logical vector of length ceil(log2(numSec)) that marks
                % which powers of 2 to keep in order to reconstruct the whole ring
                flags = split(dec2bin(numSec),'');
                flags = flags(~cellfun(@isempty,flags));
                flags = flipud(logical(cellfun(@str2num, flags)));

                tls = tls0{cur_lvl};
                if isempty(tls) % means we have to update the element for this R
                    tls = gen_gds_structure(opts, cur_lvl, R);
                    tls = tls{1};
                    tls.sname = ['EL_', num2str(cur_lvl), '_R_', num2str(ir)];
                    allels{ir} = tls;
                    trans = [0 0];
                else
                    trans = [R 0];
                end

                % initial unit cell section
                ring = cell(1, sum(2*flags)-flags(1));
                iring = 1;
                st.angle = 0;
                init = gds_element('sref', 'sname', sname(tls), 'xy',trans, 'strans',st);
                if flags(1)
                    st.angle = alpha/2;
                    ring{iring} = gds_element('sref', 'sname', sname(tls), 'xy',rotd(trans, alpha/2), 'strans',st);
                    iring = iring + 1;
                end
                sec = gds_structure(['SEC_0_R_', num2str(ir)], init);
                sects = cell(1, length(flags));
                sects{1} = sec;
                for lvl = 1:length(flags)-1
                    st.angle = alpha * 2.^(lvl-2);
                    copy1 = gds_element('sref', 'sname', sname(sec), 'xy',[0 0], 'strans',st);
                    st.angle = -st.angle;
                    copy2 = gds_element('sref', 'sname', sname(sec), 'xy',[0 0], 'strans',st); 
                    sec = gds_structure(['SEC_', num2str(lvl), '_R_', num2str(ir)], {copy1, copy2}); 
                    sects{lvl+1} = sec;
                    if flags(lvl+1)
                        st.angle = alpha * (2.^(lvl-1)+sum(pow2(find(flags(1:lvl))-1)));
                        ring{iring} = gds_element('sref', 'sname', sname(sec), 'xy',[0 0], 'strans', st);
                        iring = iring + 1;
                    end
                end
            areaFraction = (R/opts.Rmax)^2;
            allsects{ir} = sects; 
            ir = ir+1; % incriment before the assignment, becasue of the center hole
            rings{ir} = ring; 
            if areaFraction - areaFraction_old >= 0.01
                waitbar(areaFraction, wb, sprintf('%.0f%% done', areaFraction*100));
                areaFraction_old = areaFraction;
            end
            R = R + PR;
            end % while loop
            % flatten the cell arrays
            rings = [rings{:}];
            allsects = [allsects{:}];
        case 'cyl'
            %%% algorithm for rectangular patterning (cylindrical lens, deflector)
            %%% simply use the set of stripes and translate them in X- (R-)
            %%% direction. These stripes are kept in allsects cell array.
            R = opts.Rmin; % R is the inner radius of the current ring
            while R < opts.Rmax(1)
                % generating values for phase each time, in case the ring sizes are unequal
                cur_lvl = get_el_ind(opts, phase_data(R+PR/2)); % see below for definition
                PR = opts.period(cur_lvl, 1);
                tls = allsects{cur_lvl};
                rings{ir} = gds_element('sref', 'sname', sname(tls), 'xy',[R 0], 'strans',struct('angle',0));
                areaFraction = R/opts.Rmax(1);
                if areaFraction - areaFraction_old >= 0.01
                    waitbar(areaFraction, wb, sprintf('%.0f%% done', areaFraction*100));
                    areaFraction_old = areaFraction;
                end
                ir = ir + 1;
                R = R + PR;
            end
    end % switch type
    % remove empty cells
    rings(cellfun(@isempty, rings)) = [];
    allsects(cellfun(@isempty, allsects)) = [];
    allels(cellfun(@isempty, allels)) = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Make the label
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isempty(opts.label.string)
        waitbar(1, wb, 'Making the label...');
        a = opts.label.max_feature/sqrt(2);
        h = opts.label.height;
        dgrid = 2*a;
        txt = gdsii_boundarytext(opts.label.string, [0,0], h, opts.label.angle, 6);
        stxt = struct(txt);
        py = round(h/dgrid);
        txt_len = max(stxt.data.xy{end}(:,1));
        px = round(txt_len/dgrid);
        % square grid
        [X,Y] = ndgrid(0:px,0:py);
        % adding the second grid with an offset to make a checkerboard
        X = [X(:); reshape(X(1:end-1,1:end-1)+1/2,[],1)];
        Y = [Y(:); reshape(Y(1:end-1,1:end-1)+1/2,[],1)];
        xy = [0,0; a,0; a,a; 0,a] + reshape([X(:), Y(:)].', 1, 2, [])*dgrid;
        xy = squeeze(mat2cell(xy,size(xy,1),2,ones(1,size(xy,3))));
        be = gds_element('boundary', 'xy',xy, 'layer',opts.layer);
        txt = poly_bool(txt, be, 'and', 'layer', opts.layer);
        txt_struct = gds_structure('TXT', txt);
        
        if isnumeric(opts.label.location)
            position = opts.label.location;
        else
            switch opts.label.location
                case 'ne'
                    if contains(opts.type, 'cyl')
                        % square geometry with zero in the bottom left corner
                        position = opts.Rmax + [h -2*h];
                    else % circular, with origin in the middle
                        position = [opts.Rmax*cos(asin((opts.Rmax-min(12*h,opts.Rmax*.8))/opts.Rmax)),...
                            max(opts.Rmax-2*h, opts.Rmax*.8)];
                    end
                case 'nw'
                    if contains(opts.type, 'cyl')
                        % square geometry with zero in the bottom left corner
                        position = [-txt_len - h, opts.Rmax(2)-2*h];
                    else % circular, with origin in the middle
                        position = [-opts.Rmax*cos(asin((opts.Rmax-min(12*h,opts.Rmax*.8))/opts.Rmax)),...
                            max(opts.Rmax-2*h, opts.Rmax*.8)];
                    end
                case 'out'
                    if contains(opts.type, 'cyl')
                        position = opts.Rmax - [txt_len -h];
                    else
                        position = opts.Rmax - [txt_len h];
                    end
                case 'center'
                    if contains(opts.type, 'cyl')
                        position = opts.Rmax.*[.5 1] + [-txt_len/2 h];
                    else
                        position = [-txt_len/2, opts.Rmax + 5*h];
                    end
                otherwise
                    error('Wrong label location option.');
            end
        end
        label_ref = {gds_element('sref', 'sname', sname(txt_struct), 'xy',position)};
        if opts.label.mirror
            if contains(opts.type, 'cyl')
                position = -position + opts.Rmax;
            else
                position = -position;
            end
            label_ref = [label_ref, {gds_element('sref', 'sname', sname(txt_struct), ...
                'xy',position, 'strans',struct('angle', 180))}];
        end
        txt_struct = {txt_struct};
    else
        [label_ref, txt_struct] = deal({});
    end
    
    %%% gather everything on top and create the library
    top = gds_structure('TOP', [rings, label_ref]); 
    L = gds_library('META.DB', 'uunit',opts.uunit, 'dbunit',opts.dbunit, ...
        [allsects, txt_struct, allels, {top}]);

    %%% write on disk
    fname = opts.fname;
    if ~isempty(fname)
        % add file extenison if needed
        if length(fname) < 4 || ~strcmp(fname(end-3:end), '.gds')
            fname = [fname, '.gds'];
        end
        
        % format the name for the file for waitbar
        if fname(1) == '!'
            name = fname(2:end);
        else
            name = fname;
        end
        check = name == '_';
        if any(check)
            check = find(check);
            for ic = 1:length(check)
                name = [name(1:check(ic)-1), '\', name(check(ic):end)];
                check(ic:end) = check(ic:end) + 1;
            end
        end
        waitbar(1, wb, ['Writing ', name, ' on disk...']);
        
        % actually write the file
        write_gds_library(L, fname);
    end
    close(wb) % close the waitbar
    toc
end

function cur_lvl = get_el_ind(opts, cur_phase)
% Subroutine to get the index of the unit cell to use for
% the current phase value. Written as its own function 
% to separate code for 'cyl' and 'circ'
    cur_phase = mod(cur_phase, 2*pi);
    [~,cur_lvl] = min(abs(opts.phi0-cur_phase));
    % not a particularly pretty way to wrap the phase around 2*pi, e.g.
    % phase of 0.99 (units 2*pi) is closer to opts.phi0=0.01 than 0.95
    if cur_lvl == 1 && abs(opts.phi0(end)-2*pi-cur_phase) < abs(opts.phi0(1)-cur_phase)
        cur_lvl = numel(opts.phi0);
    elseif cur_lvl == numel(opts.phi0) && abs(opts.phi0(1)+2*pi-cur_phase) < abs(opts.phi0(end)-cur_phase)
        cur_lvl = 1;
    end
end