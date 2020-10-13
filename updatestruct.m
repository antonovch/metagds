function [sold, rem] = updatestruct(sold, snew, newfieldopt)
% updates old struct (sold) using fields from the new one (snew) 
% snew can be given as a cell array of name/value pairs
% rem - remainder struct with fields from snew that are not in sold. If
% snew was a cell array, so becomes rem.
% newfieldopt - option that says what to do with fields in snew, that are not in sold 
%   'ignore' -- do not copy to sold, but save in rem
%   'copy' -- copy to sold
%   'error' -- throw an error (typos, bugs, etc)
% if nargout == 1, the default is 'error', if 2 -- 'ignore'
% NOTES: not implemented for struct arrays

    if nargin < 3
        if nargout == 1
        	newfieldopt = 'error';
        elseif nargout == 2
            newfieldopt = 'ignore';
        end
    end
    was_cell = false;
    if iscell(snew)
        if isempty(snew)
            rem = {};
            return
        elseif isstruct(snew{1}) || length(snew) == 1
            snew = snew{1};
        else
            % first make sure all the field names are unique (use the last
            % value provided for diplicates) 
            names = snew(1:2:end);
            vals = snew(2:2:end);
            [names, ~, ic] = unique(names, 'stable');
            valsnew = cell(1,length(names));
            for ii = 1:length(names)
                f = find(ic == ii);
                valsnew(ii) = vals(f(end));
            end
            snew = [names;valsnew];
            snew = snew(:).';
            snew = struct(snew{:});
            was_cell = true;
        end
    end
    rem = struct();
    fnames = fieldnames(snew);
    for ii = 1:length(fnames)
        if isfield(sold, fnames{ii})
            if isstruct(sold.(fnames{ii})) && isstruct(snew.(fnames{ii}))
                sold.(fnames{ii}) = updatestruct(sold.(fnames{ii}), ...
                    snew.(fnames{ii}), newfieldopt);
            else
                sold.(fnames{ii}) = snew.(fnames{ii});
            end
        else
            switch newfieldopt
                case 'copy'
                    sold.(fnames{ii}) = snew.(fnames{ii});
                case 'ignore'
                    rem.(fnames{ii}) = snew.(fnames{ii});
                case 'error'
                    error('Wrong option %s', fnames{ii});
            end
        end
    end
    if was_cell
        rem = struct2nv(rem);
    end
end