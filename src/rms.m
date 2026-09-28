function y = rms(x, varargin)

if nargin >= 2 && ~isempty(varargin{1})
    dim = varargin{1};
else
    dim = find(size(x) ~= 1, 1);
    if isempty(dim)
        dim = 1;
    end
end

y = sqrt(mean(abs(x).^2, dim));
end
