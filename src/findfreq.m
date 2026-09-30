function freq = findfreq(sig,fs)

    if(nargin < 2)
        fs = 1;
    end

    if fs <= 0
        error('findfreq:invalidFs', 'Sampling frequency fs must be positive.');
    end

    if(~isreal(sig) || ~isreal(fs))
        error('findfreq:invalidInput', 'Inputs must be real numbers.');
    end

    [~,freq,~,~,~] = sinfit(sig);

    freq = freq*fs;

end
