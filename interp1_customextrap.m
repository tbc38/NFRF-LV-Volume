function interp_vec = interp1_customextrap(varargin)
% Function to have different constant extrapolated values at the beginning
% vs the end of the interpolated vector
% will just hand this to the interp1d function, making extrapolated values
% NaN.
% Then, will fill in the extrapolated values at the beginning with the first extrap value, 
% and at the extrapolated values at the end, with the second extrap value,
% Return the interpolated vector with custom constant bookends extrap
% values

if nargin < 3
    error('Not enough input arguments. At least x, v, and xi are required.');
end

% Remove the last argument
newArgs = varargin(1:end-1);

% Call interp1 with the modified arguments list
[varargout{1:nargout}] = interp1(newArgs{:});

interp_vec = varargout{1};

extrap_val = varargin{end};

if isempty(extrap_val)
    return
elseif isscalar(extrap_val)
    idxNaN = isnan(interp_vec);
    interp_vec(idxNaN) = extrap_val;
    return;
elseif numel(extrap_val) == 2
    % replace NaNs in first half of vec
    idxNaN = isnan(interp_vec);
    idxNaN(round(length(idxNaN)/2):end) = false;
    interp_vec(idxNaN) = extrap_val(1);

    % replace NaNs in second half of vec
    idxNaN = isnan(interp_vec);
    interp_vec(idxNaN) = extrap_val(2);
    return;
end

error('Unsupported number of extrapvals')
end
