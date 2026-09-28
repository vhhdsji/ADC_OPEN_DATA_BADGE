function b = findbin(fs, fin, n)

    if nargin < 3
        error('findbin:notEnoughInputs', ...
              'Three input arguments required: fs, fin, n.');
    end

    if ~isvector(fin)
        error('findbin:invalidInput', ...
              'fin must be a vector.');
    end

    if ~isscalar(fs) || ~isscalar(n)
        error('findbin:invalidInput', ...
              'fs and n must be scalars.');
    end

    if fs <= 0 
        error('findbin:invalidFrequency', ...
              'fs must be positive.');
    end

    if n <= 0 || floor(n) ~= n
        error('findbin:invalidN', ...
              'FFT size n must be a positive integer.');
    end

    if ~isreal(fs) || ~isreal(fin) || ~isreal(n)
        error('findbin:invalidInput', ...
              'All inputs must be real numbers.');
    end

    b = zeros(size(fin));

    for i = 1:length(fin)
        bin_start = floor(fin(i) / fs * n);

        d = 0;
        while true
            b_upper = bin_start + d;
            if gcd(b_upper, n) == 1
                b(i) = b_upper;
                break;
            end

            b_lower = bin_start - d;
            if b_lower > 0 && gcd(b_lower, n) == 1
                b(i) = b_lower;
                break;
            end

            d = d + 1;
        end
    end

end
