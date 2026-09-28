function [fitout,freq,mag,dc,phi] = sinfit(sig,varargin)

    if ~isnumeric(sig) || ~isreal(sig)
        error('sinfit:invalidInput', 'Input signal must be a real numeric array.');
    end

    if isempty(sig)
        error('sinfit:emptyInput', 'Input signal cannot be empty.');
    end

    [N,M] = size(sig);
    if(N == 1)
       sig = sig';
       N = M;
    end
    sig = mean(sig,2);

    p = inputParser;
    addOptional(p, 'f0', 0, @(x) isnumeric(x) && isscalar(x) && (x >= 0) && (x <= 0.5));
    addOptional(p, 'tol', 1e-12, @(x) isnumeric(x) && isscalar(x) && (x > 0));
    addOptional(p, 'rate', 0.5, @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x <= 1));
    addOptional(p, 'fsearch', 0, @(x) isnumeric(x) && isscalar(x));
    addOptional(p, 'verbose', 0, @(x) isnumeric(x) && isscalar(x) && ismember(x, [0, 1]));
    addOptional(p, 'niter', 100, @(x) isnumeric(x) && isscalar(x) && (x > 0));
    parse(p, varargin{:});
    f0 = p.Results.f0;
    tol = p.Results.tol;
    rate = p.Results.rate;
    fsearch = p.Results.fsearch;
    verbose = p.Results.verbose;
    niter = p.Results.niter;

    if(f0 == 0)
        fsearch = 1;
        spec = abs(fft(sig));
        spec(1) = 0;
        spec = spec(1:floor(N/2));

        [~,k0] = max(spec);

        nspec = floor(N/2);
        if(spec(min(max(k0+1,1),nspec)) > spec(min(max(k0-1,1),nspec)))
            r = 1;
        else
            r = -1;
        end

        k_neighbor = min(max(k0+r, 1), nspec);
        f0 = (k0-1 + r*spec(k_neighbor)/(spec(k0)+spec(k_neighbor)))/N;
    end

    time = (0:N-1)';
    theta = 2*pi*f0*time;
    M = [cos(theta), sin(theta), ones([N,1])];
    x = linsolve(M,sig);
    A = x(1);
    B = x(2);
    dc = x(3);

    freq = f0;

    if(fsearch)
        delta_f = 0;

        for ii = 1:niter
            freq = freq+delta_f;
            theta = 2*pi*freq*time;

            M = [cos(theta), sin(theta), ones([N,1]), (-A*2*pi*time.*sin(theta)+B*2*pi*time.*cos(theta))/N];
            x = linsolve(M,sig);
            A = x(1);
            B = x(2);
            dc = x(3);
            delta_f = x(4)*rate/N;

            relerr = rms(x(4)/N*M(:,4)) / sqrt(x(1)^2+x(2)^2);

            if verbose
                fprintf('Freq iterating (%d): freq = %d, delta_f = %d, rel_err = %d\n', ii, freq, delta_f, relerr);
            end

            if(relerr < tol)
                break;
            end
        end

        if ii == niter && relerr >= tol
            warning('sinfit:noConvergence', ...
                'Failed to converge in %d iterations. Relative error = %.2e', niter, relerr);
        end
    end

    fitout = A*cos(theta)+B*sin(theta)+dc;

    mag = sqrt(A^2+B^2);
    phi = -atan2(B,A);
end
