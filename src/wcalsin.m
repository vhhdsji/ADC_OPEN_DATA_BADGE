function [weight,offset,postcal,ideal,err,freqcal] = wcalsin(bits,varargin)

    if iscell(bits)

        ND = numel(bits);
        if ND == 0
            error('wcalsin:EmptyInput','Empty cell array for bits.');
        end

        p = inputParser;
        addOptional(p, 'freq', 0, @(x) isnumeric(x) && isvector(x) && all(x>=0) && all(x<=0.5));
        addOptional(p, 'rate', 0.5, @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x < 1));
        addOptional(p, 'reltol', 1E-12, @(x) isnumeric(x) && isscalar(x) && (x > 0));
        addOptional(p, 'niter', 100, @(x) isnumeric(x) && isscalar(x) && (x > 0));
        addOptional(p, 'order', 1, @(x) isnumeric(x) && isscalar(x) && (x > 0));
        addOptional(p, 'fsearch', 0, @(x) isnumeric(x) && isscalar(x));
        addOptional(p, 'verbose', 0, @(x) (islogical(x) && isscalar(x)) || (isnumeric(x) && isscalar(x) && ismember(x, [0, 1])));
        addOptional(p, 'autotrans', 1, @(x) isnumeric(x) && isscalar(x) && ismember(x, [0, 1]));
        addOptional(p, 'autopatch', 1, @(x) isnumeric(x) && isscalar(x) && ismember(x, [0, 1]));
        addOptional(p, 'nomWeight', []);
        parse(p, varargin{:});
        freq = p.Results.freq;
        order = max(round(p.Results.order),1);
        nomWeight = p.Results.nomWeight;
        rate = p.Results.rate;
        reltol = p.Results.reltol;
        niter = p.Results.niter;
        fsearch = p.Results.fsearch;
        verbose = p.Results.verbose;
        autotrans = p.Results.autotrans;
        autopatch = p.Results.autopatch;

        bits_cell = cell(1,ND);
        Nk = zeros(1,ND);
        Mk = zeros(1,ND);
        for k = 1:ND
            Bk = bits{k};
            if isempty(Bk)
                error('wcalsin:EmptyDataset','Dataset %d is empty.',k);
            end
            [nTmp,mTmp] = size(Bk);
            if autotrans && nTmp < mTmp
                Bk = Bk';
                [nTmp,mTmp] = size(Bk);
            end
            bits_cell{k} = Bk;
            Nk(k) = nTmp;
            Mk(k) = mTmp;
        end
        if any(Mk ~= Mk(1))
            error('wcalsin:InconsistentWidth','All datasets must have the same number of columns (bits).');
        end
        M_orig = Mk(1);
        if isempty(nomWeight)
            nomWeight = 2.^(M_orig-1:-1:0);
        end

        if isscalar(freq)
            freq = ones(1,ND) * freq;
        end
        if numel(freq) ~= ND
            error('wcalsin:FreqLength','Length of freq vector must match number of datasets.');
        end

        for k = 1:ND
            if freq(k) == 0 || fsearch == 1
                [~,~,~,~,~,fk] = wcalsin(bits_cell{k},freq(k),rate,reltol,niter,order,1,verbose,0,autopatch,nomWeight);
                freq(k) = fk;
            end
        end

        bits_all = vertcat(bits_cell{:});
        Ntot = size(bits_all,1);

        Lmap = 1:M_orig;
        Kmap = ones(1,M_orig);

        if autopatch && rank([bits_all,ones(Ntot,1)]) < M_orig+1
            warning('Rank deficiency detected across datasets. Try patching...');
            bits_patch_all = [];
            LR = [];
            M2 = 0;
            for i1 = 1:M_orig
                if max(bits_all(:,i1))==min(bits_all(:,i1))
                    Lmap(i1) = 0;
                    Kmap(i1) = 0;
                    if verbose
                        fprintf('Constant column %d discarded\n', i1);
                    end
                elseif(rank([ones(Ntot,1),bits_patch_all,bits_all(:,i1)]) > rank([ones(Ntot,1),bits_patch_all]))
                    bits_patch_all = [bits_patch_all,bits_all(:,i1)];
                    LR = [LR,i1];
                    [~,M2] = size(bits_patch_all);
                    Lmap(i1) = M2;
                else
                    flag = 0;
                    for i2 = 1:M2
                        r1 = bits_all(:,i1)-mean(bits_all(:,i1));
                        r2 = bits_patch_all(:,i2)-mean(bits_patch_all(:,i2));
                        cor = mean(r1.*r2)/local_rms(r1)/local_rms(r2);
                        if(abs(abs(cor)-1) < 1E-3)
                            Lmap(i1) = i2;
                            Kmap(i1) = nomWeight(i1)/nomWeight(LR(i2));
                            bits_patch_all(:,i2) = bits_patch_all(:,i2) + bits_all(:,i1)*Kmap(i1);
                            if verbose
                                fprintf('Patched column %d -> column %d (ratio: %.6g)\n', i1, LR(i2), Kmap(i1));
                            end
                            flag = 1;
                            break;
                        end
                    end
                    if(flag == 0)
                        Lmap(i1) = 0;
                        Kmap(i1) = 0;
                        warning('Patch warning: cannot find correlated column for column %d. Resulting weight will be zero',i1);
                    end
                end
            end
            [~,M_patch] = size(bits_patch_all);
            if(rank([ones(Ntot,1),bits_patch_all]) < M_patch+1)
                error('Patch failed: rank still deficient after patching across datasets. Try adjusting nomWeight.');
            end
        else
            bits_patch_all = bits_all;
            M_patch = M_orig;
        end

        if M_patch == 0
            error('Patched bits are empty. No valid columns remain after patching.');
        end

        MAG = floor(log10(max(abs([max(bits_patch_all);min(bits_patch_all)]))));
        MAG(isinf(MAG)) = 0;
        bits_patch_all = bits_patch_all.*(ones(Ntot,1)*10.^(-MAG));

        numHCols = ND*order;
        xc = zeros(Ntot, numHCols);
        xs = zeros(Ntot, numHCols);

        rowStart = 1;
        for k = 1:ND
            Nk_k = Nk(k);
            rowEnd = rowStart + Nk_k - 1;

            theta_mat = (0:(Nk_k-1))' * freq(k) * (1:order);
            xc(rowStart:rowEnd, (k-1)*order+1 : k*order) = cos(theta_mat*2*pi);
            xs(rowStart:rowEnd, (k-1)*order+1 : k*order) = sin(theta_mat*2*pi);

            rowStart = rowEnd + 1;
        end

        A1 = [bits_patch_all,ones(Ntot,1),xc(:,2:end),xs];
        b = -xc(:,1);
        x1 = linsolve(A1,b);

        A2 = [bits_patch_all,ones(Ntot,1),xs(:,2:end),xc];
        b = -xs(:,1);
        x2 = linsolve(A2,b);

        if(local_rms(A1*x1-b) < local_rms(A2*x2-b))
            x = x1;
            sel = 0;
        else
            x = x2;
            sel = 1;
        end

        w0 = sqrt(1 + x(M_patch + numHCols + 1)^2);

        wpatch = (x(1:M_patch)'/w0) .* (10.^-MAG);
        weight = wpatch(max(Lmap,1)) .* Kmap;
        offset = -x(M_patch+1)/w0;

        postcal = cell(1,ND);
        ideal = cell(1,ND);
        err = cell(1,ND);

        if(sel)
            ideal{1} = -(xs(1:Nk(1),1) + xs(1:Nk(1),2:order) * x(M_patch + 1 + 1 : M_patch + 1 + order-1) + xc(1:Nk(1),1:order) * x(M_patch + 1 + numHCols : M_patch + 1 + numHCols + order-1))'/w0;
        else
            ideal{1} = -(xc(1:Nk(1),1) + xc(1:Nk(1),2:order) * x(M_patch + 1 + 1 : M_patch + 1 + order-1) + xs(1:Nk(1),1:order) * x(M_patch + 1 + numHCols : M_patch + 1 + numHCols + order-1))'/w0;
        end

        rowStart = Nk(1)+1;
        for k = 2:ND
            rowEnd = rowStart + Nk(k) - 1;
            if(sel)
                ideal{k} = -(xs(rowStart:rowEnd,(k-1)*order+1:k*order) * x(M_patch + 1 + (k-1)*order : M_patch + 1 + k*order - 1) +...
                             xc(rowStart:rowEnd,(k-1)*order+1:k*order) * x(M_patch + 1 + numHCols + (k-1)*order : M_patch + 1 + numHCols + k*order - 1) )'/w0;
            else
                ideal{k} = -(xc(rowStart:rowEnd,(k-1)*order+1:k*order) * x(M_patch + 1 + (k-1)*order : M_patch + 1 + k*order - 1) +...
                             xs(rowStart:rowEnd,(k-1)*order+1:k*order) * x(M_patch + 1 + numHCols + (k-1)*order : M_patch + 1 + numHCols + k*order - 1) )'/w0;
            end
            rowStart = rowEnd+1;
        end

        for k = 1:ND
            postcal{k} = weight * bits_cell{k}';
            err{k} = postcal{k} - offset - ideal{k};
        end

        if sum(weight) < 0
            weight = -weight;
            offset = -offset;
            for k = 1:ND
                postcal{k} = -postcal{k};
                ideal{k}   = -ideal{k};
                err{k}     = -err{k};
            end
        end

        freqcal = freq;
        return;
    end

    p = inputParser;
    addOptional(p, 'freq', 0, @(x) isnumeric(x) && isscalar(x) && (x >= 0) && (x <= 0.5));
    addOptional(p, 'rate', 0.5, @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x < 1));
    addOptional(p, 'reltol', 1E-12, @(x) isnumeric(x) && isscalar(x) && (x > 0));
    addOptional(p, 'niter', 100, @(x) isnumeric(x) && isscalar(x) && (x > 0));
    addOptional(p, 'order', 1, @(x) isnumeric(x) && isscalar(x) && (x > 0));
    addOptional(p, 'fsearch', 0, @(x) isnumeric(x) && isscalar(x));
    addOptional(p, 'verbose', 0, @(x) (islogical(x) && isscalar(x)) || (isnumeric(x) && isscalar(x) && ismember(x, [0, 1])));
    addOptional(p, 'autotrans', 1, @(x) isnumeric(x) && isscalar(x) && ismember(x, [0, 1]));
    addOptional(p, 'autopatch', 1, @(x) isnumeric(x) && isscalar(x) && ismember(x, [0, 1]));
    addOptional(p, 'nomWeight', []);
    parse(p, varargin{:});
    freq = p.Results.freq;
    order = max(round(p.Results.order),1);
    nomWeight = p.Results.nomWeight;
    rate = p.Results.rate;
    reltol = p.Results.reltol;
    niter = p.Results.niter;
    fsearch = p.Results.fsearch;
    verbose = p.Results.verbose;
    autotrans = p.Results.autotrans;
    autopatch = p.Results.autopatch;

    [N,M] = size(bits);
    if autotrans && N < M
        bits = bits';
        [N,M] = size(bits);
    end
    if isempty(nomWeight)
        nomWeight = 2.^(M-1:-1:0);
    end

    L = [1:M];
    LR = [1:M];
    K = ones(1,M);

    if autopatch && rank([bits,ones(N,1)]) < M+1
        warning('Rank deficiency detected. Try patching...');
        bits_patch = [];
        LR = [];
        M2 = 0;
        for i1 = 1:M
            if(max(bits(:,i1))==min(bits(:,i1)))
                L(i1) = 0;
                K(i1) = 0;
                if verbose
                    fprintf('Constant column %d discarded\n', i1);
                end
            elseif(rank([ones(N,1),bits_patch,bits(:,i1)]) > rank([ones(N,1),bits_patch]))
                bits_patch = [bits_patch,bits(:,i1)];
                LR = [LR,i1];
                [~,M2] = size(bits_patch);
                L(i1) = M2;
            else
                flag = 0;
                for i2 = 1:M2
                    r1 = bits(:,i1)-mean(bits(:,i1));
                    r2 = bits_patch(:,i2) - mean(bits_patch(:,i2));
                    cor = mean(r1.*r2)/local_rms(r1)/local_rms(r2);
                    if(abs(abs(cor)-1) < 1E-3)
                        L(i1) = i2;
                        K(i1) = nomWeight(i1)/nomWeight(LR(i2));

                        bits_patch(:,i2) = bits_patch(:,i2) + bits(:,i1)*nomWeight(i1)/nomWeight(LR(i2));
                        if verbose
                            fprintf('Patched column %d -> column %d (ratio: %.6g)\n', i1, LR(i2), K(i1));
                        end
                        flag = 1;
                        break;
                    end
                end
                if(flag == 0)
                    L(i1) = 0;
                    K(i1) = 0;
                    warning('Patch warning: cannot find the correlated column for column %d. The resulting weight will be zero',i1);
                end
            end
        end
        [~,M] = size(bits_patch);
        if(rank([ones(N,1),bits_patch]) < M+1)
            error('Patch failed: rank still deficient after patching. This may be fixed by changing nomWeight.')
        end
    else
        bits_patch = bits;
    end

    if M == 0
        error('Patched bits are empty. No valid columns remain after patching.');
    end

    MAG = floor(log10(max(abs([max(bits_patch);min(bits_patch)]))));
    MAG(isinf(MAG)) = 0;
    bits_patch = bits_patch.*(ones(N,1)*10.^(-MAG));

    if(freq == 0)
        fsearch = 1;
        freq = [];
        for i1 = 1:min(M,5)
            if verbose
                fprintf('Freq coarse searching (%d/5):',i1);
            end

            freq = [freq, findfreq(bits_patch(:,1:i1)*nomWeight(LR(1:i1))',1)];
            if verbose
                fprintf(' freq = %d\n',freq(end));
            end
        end
        freq = median(freq);
    end

    theta_mat = (0:(N-1))'*freq*(1:order);
    xc = cos(theta_mat*2*pi);
    xs = sin(theta_mat*2*pi);

    A1 = [bits_patch(1:N,1:M),ones(N,1),xc(:,2:end),xs];
    b = -xc(:,1);
    x1 = linsolve(A1,b);

    A2 = [bits_patch(1:N,1:M),ones(N,1),xs(:,2:end),xc];
    b = -xs(:,1);
    x2 = linsolve(A2,b);
    if(local_rms(A1*x1-b) < local_rms(A2*x2-b))
        x = x1;
        sel = 0;
    else
        x = x2;
        sel = 1;
    end

    if(fsearch)

        warning_state = warning;
        warning_cleanup = onCleanup(@() warning(warning_state));
        warning('off','all');

        delta_f = 0;
        time_mat = (0:(N-1))'*ones([1,order]);

        for ii = 1:niter
            freq = freq+delta_f;
            theta_mat = (0:(N-1))'*freq*(1:order);

            xc = cos(theta_mat*2*pi);
            xs = sin(theta_mat*2*pi);

            order_mat = ones([N,1]) * (1:order);
            if(sel)
                KS = ones([N,1]) * [1,x(M+2:M+order)'] .* order_mat;
                KC = ones([N,1]) * x(M+1+order:M+order*2)' .* order_mat;
            else
                KC = ones([N,1]) * [1,x(M+2:M+order)'] .* order_mat;
                KS = ones([N,1]) * x(M+1+order:M+order*2)' .* order_mat;
            end

            xcd = -2*pi * KC .* time_mat .* sin(theta_mat*2*pi) / N;
            xsd =  2*pi * KS .* time_mat .* cos(theta_mat*2*pi) / N;

            A = [bits_patch(1:N,1:M),ones(N,1),xc(:,2:end),xs,sum(xcd+xsd,2)];
            b = -xc(:,1);
            x1 = linsolve(A,b);
            e1 = A*x1-b;

            A = [bits_patch(1:N,1:M),ones(N,1),xs(:,2:end),xc,sum(xcd+xsd,2)];
            b = -xs(:,1);
            x2 = linsolve(A,b);
            e2 = A*x2-b;

            if(local_rms(e1) < local_rms(e2))
                x = x1;
                sel = 0;
            else
                x = x2;
                sel = 1;
            end

            delta_f = x(end)*rate /N;
            relerr = local_rms(x(end)/N*A(:,end)) / sqrt(1+x(M+1+order)^2);

            if verbose
                fprintf('Freq fine iterating (%d): freq = %d, delta_f = %d, rel_err = %d\n',ii,freq,delta_f, relerr);
            end

            if(relerr < reltol)
                break;
            end

        end

        if ii == niter && relerr >= reltol
            warning('wcalsin:noConvergence', ...
                'Failed to converge in %d iterations. Relative error = %.2e', niter, relerr);
        end

        clear warning_cleanup
    end

    w0 = sqrt(1+x(M+1+order)^2);

    weight = x(1:M)'/w0.*(10.^-MAG);
    weight = weight(max(L,1)).*K;

    offset = -x(M+1)/w0;

    postcal = weight*bits';

    if(sel)

        ideal = -(xs(:,1) + xs(:,2:end) * x(M+2:M+order) + xc * x(M+1+order:M+2*order))'/w0;
    else

        ideal = -(xc(:,1) + xc(:,2:end) * x(M+2:M+order) + xs * x(M+1+order:M+2*order))'/w0;
    end

    err = postcal-offset-ideal;

    if(sum(weight)<0)
        weight = -weight;
        offset = -offset;
        postcal = -postcal;
        ideal = -ideal;
        err = -err;
    end

    freqcal = freq;

    snr_linear = std(ideal) / std(err);
    if snr_linear < 10
        warning('SNR (%.1f dB) is below 20 dB. Calibration may have failed or sinewave may not be correctly extracted.', 20*log10(snr_linear));
    end

end
