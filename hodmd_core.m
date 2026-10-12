function varargout = hodmd_core(action,varargin)
% Shared numerical functions matching step02_run_HODMD.m.
% The wrapper only connects the existing numerical functions without plotting.
switch action
    case 'defaults'
        cfg.pre=struct('smoothAlongQ',true,'qSmoothWindow',20);
        cfg.hodmd=struct('delayOrder',20,'dmdRank',9);
        cfg.target=struct('period_fs',303.24,'tolerance_fs',20, ...
            'numPositiveModesPerTarget',1);
        cfg.notes='I2 preprocessing: q moving average only; no coordinate reduction.';
        varargout={cfg};
    case 'formFactor'
        varargout={i2_atomic_form_factor(varargin{1})};
    case 'spectrum'
        [q,t,A,lambda,ids]=deal(varargin{:});
        h=struct('lambda',lambda);
        varargout={compute_hodmd_q_spectrum(q,t,A,h,ids)};
    case {'analyze','component'}
        [A,q,t,cfg]=deal(varargin{:});
        assert(all(isfinite(A(:))) && isreal(A),'HODMD input must be finite and real.');
        assert(isequal(size(A),[numel(q),numel(t)]),'HODMD input dimensions disagree.');
        dt=median(diff(t)); d=cfg.hodmd.delayOrder;
        assert(dt>0 && max(abs(diff(t)-dt))<1e-6*dt,'HODMD requires uniform time.');
        assert(d==0 || (d>=2 && d<=numel(t)-2),'Invalid HODMD delay order.');
        B=preprocess_hodmd_matrix(A,q,cfg);
        Uq=eye(size(A,1));
        Z=make_delay_matrix(B,d);
        [h,di]=exact_dmd(Z(:,1:end-1),Z(:,2:end),cfg.hodmd.dmdRank,dt);
        [h.beta,h.b_firstSnapshot,h.betaFitInfo]=fit_global_dmd_amplitudes(h.Phi,h.lambda,Z);
        h.b=h.beta;
        h=sort_hodmd_modes_by_amplitude(h);
        mt=build_hodmd_mode_table(h,dt);
        groups=find_target_hodmd_modes(h,mt,cfg.target);
        ids=unique([groups{:}],'stable');
        rec=Uq*reconstruct_hodmd_modes(h,ids,size(B,1),d,numel(t));
        if strcmp(action,'component')
            varargout={rec,struct('selectedIndices',ids,'modeTable',mt)};
            varargout=varargout(1:nargout);
            return;
        end
        spec=compute_hodmd_q_spectrum(q,t,rec,h,ids);
        out=struct('modeTable',mt,'selectedIndices',ids,'reconstruction',rec, ...
            'spectrum',spec,'lambda',h.lambda,'beta',h.beta, ...
            'dmdInfo',di,'betaFitInfo',h.betaFitInfo,'cfg',cfg);
        candidates=find(mt.Frequency>0 & abs(mt.Period_fs-cfg.target.period_fs)<=cfg.target.tolerance_fs);
        out.nearbyModeTable=mt(candidates,:);
        out.nearbyPairCount=numel(candidates);
        out.nearbyReconstructions=cell(numel(candidates),1);
        for j=1:numel(candidates)
            ii=include_conjugate_pairs(h.lambda,candidates(j));
            y=Uq*reconstruct_hodmd_modes(h,ii,size(B,1),d,numel(t));
            out.nearbyReconstructions{j}=y;
        end
        varargout={out};
    case 'audit'
        original=fileread(varargin{1}); own=fileread([mfilename('fullpath') '.m']);
        names={'preprocess_hodmd_matrix','make_delay_matrix','fold_delay_matrix_to_series','exact_dmd','fit_global_dmd_amplitudes','sort_hodmd_modes_by_amplitude','apply_hodmd_sort_index','build_hodmd_mode_table','find_target_hodmd_modes','include_conjugate_pairs','find_conjugate_partner','unique_valid_indices','reconstruct_hodmd_modes','compute_hodmd_q_spectrum','build_modal_time_basis','append_oscillatory_basis','make_basis_info','i2_atomic_form_factor'};
        matched=false(numel(names),1);
        for j=1:numel(names)
            matched(j)=strcmp(extract_function(original,names{j}),extract_function(own,names{j}));
        end
        varargout={table(string(names(:)),matched,'VariableNames',{'Function','MatchesOriginal'})};
    otherwise
        error('Unknown adapter action: %s',action);
end
end

function text=extract_function(source,name)
    starts=regexp(source,'(?m)^function\>');
    starts=[starts,numel(source)+1]; text='';
    for k=1:numel(starts)-1
        part=source(starts(k):starts(k+1)-1);
        head=part(1:find(part=='(',1)-1);
        if ~isempty(regexp(head,['\<' name '\s*$'],'once'))
            text=strtrim(strrep(part,sprintf('\r\n'),sprintf('\n'))); return;
        end
    end
    assert(~isempty(text),'Missing source function: %s',name);
end
function A_work = preprocess_hodmd_matrix(A_in, q, cfg)
    assert(size(A_in, 1) == numel(q), 'q size does not match A.');
    assert(all(isfinite(A_in(:))), 'A must contain finite values.');
    A_work = A_in;
    if cfg.pre.smoothAlongQ
        w = cfg.pre.qSmoothWindow;
        assert(isscalar(w) && isfinite(w) && w >= 1 && w == floor(w), ...
            'qSmoothWindow must be a positive integer.');
        A_work = movmean(A_work, min(w, size(A_work, 1)), 1, 'Endpoints', 'shrink');
    end
end


function Z = make_delay_matrix(Y, d)
    % Y is the HODMD/DMD input coordinate matrix B: r x Nt
    % If d = 0, no delay embedding is applied and Z = Y.
    %
    % z_k = [Y(:,k);
    %        Y(:,k+1);
    %        ...
    %        Y(:,k+d-1)]
    %
    % Z: r*d x K, K = Nt-d+1

    [r, Nt] = size(Y);

    if d == 0
        Z = Y;
        return;
    end

    K = Nt - d + 1;

    Z = zeros(r*d, K);

    for k = 1:K
        block = Y(:, k:k+d-1);
        Z(:, k) = block(:);
    end
end

function Y = fold_delay_matrix_to_series(Z, r, d, Nt)
    % Inverse of delay embedding by diagonal averaging.
    %
    % Z: r*d x K
    % Y: r x Nt

    K = size(Z, 2);

    Y = zeros(r, Nt);
    W = zeros(1, Nt);

    for k = 1:K
        z = Z(:, k);
        block = reshape(z, r, d);

        for ell = 1:d
            tidx = k + ell - 1;
            Y(:, tidx) = Y(:, tidx) + block(:, ell);
            W(tidx) = W(tidx) + 1;
        end
    end

    Y = Y ./ W;
end

function [hodmd, info] = exact_dmd(Z1, Z2, rDMD, dt)
    % Exact DMD:
    %   Z2 ≈ A Z1
    %
    % Returns eigenvalues lambda and exact DMD modes Phi.

    [U, S, V] = svd(Z1, 'econ');
    s = diag(S);

    r = min([rDMD, numel(s), size(U,2)]);

    U_r = U(:, 1:r);
    S_r = S(1:r, 1:r);
    V_r = V(:, 1:r);

    Atilde = U_r' * Z2 * V_r / S_r;

    [W, D] = eig(Atilde);

    lambda = diag(D);

    Phi = Z2 * V_r / S_r * W;

    hodmd = struct();
    hodmd.lambda = lambda;
    hodmd.Phi = Phi;
    hodmd.W = W;
    hodmd.Atilde = Atilde;
    hodmd.dt = dt;

    info = struct();
    info.delaySVDValues = s;
    info.delayCumulativeEnergy = cumsum(s.^2) / sum(s.^2);
    info.dmdRank = r;
end

function [beta, bFirst, fitInfo] = fit_global_dmd_amplitudes(Phi, lambda, Z)
    % Fit beta from all delay/direct snapshots:
    %   Z(:,k) ~= sum_j Phi(:,j) * beta_j * lambda_j^(k-1)
    %
    % This is a global least-squares amplitude fit over the full trajectory.

    lambda = lambda(:);
    nModes = numel(lambda);
    K = size(Z, 2);

    timePowers = zeros(nModes, K);

    for j = 1:nModes
        timePowers(j, :) = lambda(j).^(0:K-1);
    end

    X = zeros(numel(Z), nModes);

    for j = 1:nModes
        X(:, j) = reshape(Phi(:, j) * timePowers(j, :), [], 1);
    end

    zVec = Z(:);
    beta = X \ zVec;
    bFirst = Phi \ Z(:, 1);

    zFitGlobal = X * beta;
    zFitFirst = X * bFirst;

    sX = svd(X, 'econ');

    if isempty(sX) || sX(1) == 0
        condX = inf;
        rankX = 0;
    else
        tolX = max(size(X)) * eps(sX(1));
        rankX = sum(sX > tolX);

        if sX(end) <= tolX
            condX = inf;
        else
            condX = sX(1) / sX(end);
        end
    end

    zNorm = max(norm(zVec), eps);

    fitInfo = struct();
    fitInfo.method = 'global_delay_vector_least_squares';
    fitInfo.relErr_globalBeta = norm(zVec - zFitGlobal) / zNorm;
    fitInfo.relErr_firstSnapshotB = norm(zVec - zFitFirst) / zNorm;
    fitInfo.condX = condX;
    fitInfo.rankX = rankX;
    fitInfo.numModes = nModes;
    fitInfo.numSnapshots = K;
    fitInfo.numEquations = numel(Z);
    fitInfo.betaNorm = norm(beta);
    fitInfo.bFirstNorm = norm(bFirst);
    fitInfo.betaMinusFirstNorm = norm(beta - bFirst);
    fitInfo.timePowers = timePowers;
end

function [hodmd, idx] = sort_hodmd_modes_by_amplitude(hodmd)
    b = hodmd.b;
    Phi = hodmd.Phi;

    amp = abs(b(:)) .* sqrt(sum(abs(Phi).^2, 1)).';

    [~, idx] = sort(amp, 'descend');
    hodmd = apply_hodmd_sort_index(hodmd, idx);
end

function hodmd = apply_hodmd_sort_index(hodmd, idx)
    hodmd.lambda = hodmd.lambda(idx);
    hodmd.Phi = hodmd.Phi(:, idx);

    if isfield(hodmd, 'b')
        hodmd.b = hodmd.b(idx);
    end

    if isfield(hodmd, 'beta')
        hodmd.beta = hodmd.beta(idx);
    end

    if isfield(hodmd, 'b_firstSnapshot')
        hodmd.b_firstSnapshot = hodmd.b_firstSnapshot(idx);
    end

    if isfield(hodmd, 'W')
        hodmd.W = hodmd.W(:, idx);
    end
end

function T = build_hodmd_mode_table(hodmd, dt_fs)
    lambda = hodmd.lambda(:);
    b = hodmd.b(:);
    Phi = hodmd.Phi;

    nModes = numel(lambda);

    lambdaSafe = lambda;
    small = abs(lambdaSafe) < eps;
    lambdaSafe(small) = eps;

    omega = log(lambdaSafe) / dt_fs;

    growthRate = real(omega);
    angularFreq = imag(omega);
    freq = angularFreq / (2*pi);

    period_fs = inf(nModes,1);
    validFreq = abs(freq) > eps;
    period_fs(validFreq) = 1 ./ abs(freq(validFreq));

    % Convert DMD frequency to spectroscopic wavenumber [cm^-1].
    % freq is in cycles/fs.
    c_cm_per_fs = 2.99792458e-5;
    wavenumber_cm1 = abs(freq) / c_cm_per_fs;

    decayTime_fs = nan(nModes,1);

    stability = cell(nModes,1);

    for i = 1:nModes
        if growthRate(i) < -1e-12
            decayTime_fs(i) = -1 / growthRate(i);
            stability{i} = 'decay';
        elseif growthRate(i) > 1e-12
            decayTime_fs(i) = -1 / growthRate(i);
            stability{i} = 'growth';
        else
            decayTime_fs(i) = inf;
            stability{i} = 'neutral';
        end
    end

    amp = abs(b) .* sqrt(sum(abs(Phi).^2, 1)).';

    Mode = (1:nModes).';
    AbsLambda = abs(lambda);
    AngleLambda = angle(lambda);
    Frequency = freq;
    Period_fs = period_fs;
    Wavenumber_cm1 = wavenumber_cm1;
    GrowthRate = growthRate;
    DecayTime_fs = decayTime_fs;
    Amplitude = amp;
    Stability = stability;

    T = table(Mode, AbsLambda, AngleLambda, Frequency, Period_fs, ...
        Wavenumber_cm1, GrowthRate, DecayTime_fs, Amplitude, Stability);
end

function groups = find_target_hodmd_modes(hodmd, modeTable, targetCfg)
    % Find positive-frequency modes near target periods,
    % then include conjugate pair.

    lambda = hodmd.lambda;

    nTargets = numel(targetCfg.period_fs);
    groups = cell(nTargets, 1);

    for i = 1:nTargets
        target = targetCfg.period_fs(i);
        tol = targetCfg.tolerance_fs;

        isCandidate = ...
            modeTable.Frequency > 0 & ...
            isfinite(modeTable.Period_fs) & ...
            abs(modeTable.Period_fs - target) <= tol;

        idx = find(isCandidate);

        if isempty(idx)
            groups{i} = [];
            continue;
        end

        % Sort candidates by amplitude
        [~, order] = sort(modeTable.Amplitude(idx), 'descend');
        idx = idx(order);

        nPick = min(targetCfg.numPositiveModesPerTarget, numel(idx));
        idx = idx(1:nPick);

        idxWithPairs = include_conjugate_pairs(lambda, idx);

        groups{i} = unique_valid_indices(idxWithPairs, numel(lambda));
    end
end

function idxOut = include_conjugate_pairs(lambda, idxIn)
    idxIn = idxIn(:).';
    idxOut = idxIn;

    n = numel(lambda);
    tol = 1e-6;

    for idx = idxIn
        if idx < 1 || idx > n
            continue;
        end

        pairIdx = find_conjugate_partner(lambda, idx, tol);

        if ~isempty(pairIdx)
            idxOut(end+1) = pairIdx; %#ok<AGROW>
        end
    end

    idxOut = unique(idxOut, 'stable');
end

function pairIdx = find_conjugate_partner(lambda, idx, tol)
    lambda = lambda(:);
    pairIdx = [];

    if idx < 1 || idx > numel(lambda)
        return;
    end

    lam = lambda(idx);

    if abs(imag(lam)) <= tol * max(1, abs(lam))
        return;
    end

    target = conj(lam);

    dist = abs(lambda - target);
    dist(idx) = inf;

    [minDist, candidateIdx] = min(dist);

    if isfinite(minDist) && minDist <= tol * max(1, abs(target))
        pairIdx = candidateIdx;
    end
end

function idx = unique_valid_indices(idx, nMax)
    idx = unique(idx(:).', 'stable');
    idx = idx(idx >= 1 & idx <= nMax);
end

function Y_rec = reconstruct_hodmd_modes(hodmd, modeIdx, r, d, Nt)
    % Reconstruct selected DMD/HODMD modes in reduced coordinate space.
    %
    % For d = 0, this is standard DMD reconstruction with no delay folding.
    % For d > 0, reconstruction is done in delay space and folded back
    % by diagonal averaging.

    lambda = hodmd.lambda(:);
    Phi = hodmd.Phi;
    b = hodmd.b(:);

    modeIdx = unique_valid_indices(modeIdx, numel(lambda));

    if d == 0
        Y_rec = zeros(r, Nt);

        if isempty(modeIdx)
            return;
        end

        Phi_sel = Phi(:, modeIdx);
        lambda_sel = lambda(modeIdx);
        b_sel = b(modeIdx);

        for k = 1:Nt
            timePower = lambda_sel.^(k-1);
            coef = b_sel .* timePower;
            Y_rec(:, k) = Phi_sel * coef;
        end

        Y_rec = real(Y_rec);
        return;
    end

    K = Nt - d + 1;

    Z_rec = zeros(r*d, K);

    if isempty(modeIdx)
        Y_rec = zeros(r, Nt);
        return;
    end

    Phi_sel = Phi(:, modeIdx);
    lambda_sel = lambda(modeIdx);
    b_sel = b(modeIdx);

    for k = 1:K
        timePower = lambda_sel.^(k-1);
        coef = b_sel .* timePower;
        Z_rec(:, k) = Phi_sel * coef;
    end

    Y_rec = fold_delay_matrix_to_series(Z_rec, r, d, Nt);
    Y_rec = real(Y_rec);
end

function qSpectrum = compute_hodmd_q_spectrum(q, t, A_rec, hodmd, modeIdx)
    modeIdx = unique_valid_indices(modeIdx, numel(hodmd.lambda));

    qSpectrum = struct();
    qSpectrum.valid = false;
    qSpectrum.method = 'dehankelized_reconstruction_temporal_regression';
    qSpectrum.q = q(:);
    qSpectrum.t = t(:).';
    qSpectrum.modeIdx = modeIdx;
    qSpectrum.G = [];
    qSpectrum.C = [];
    qSpectrum.basisInfo = struct([]);
    qSpectrum.condG = nan;
    qSpectrum.rankG = nan;
    qSpectrum.qCombined = [];
    qSpectrum.qModes = [];
    qSpectrum.qModeReal = [];
    qSpectrum.qModeImag = [];
    qSpectrum.qModeAbs = [];
    qSpectrum.C_real = [];
    qSpectrum.C_cos = [];
    qSpectrum.C_sin = [];
    qSpectrum.C_amp = [];
    qSpectrum.C_phase = [];
    qSpectrum.A_fit = [];
    qSpectrum.fitResidual = [];
    qSpectrum.fitResidualNorm = nan;

    if isempty(modeIdx)
        return;
    end

    q = q(:);
    t = t(:).';
    [nQ, Nt] = size(A_rec);

    if nQ ~= numel(q) || Nt ~= numel(t)
        error('A_rec must be Nq x Nt and match q and t for q-spectrum regression.');
    end

    dt = median(diff(t));
    tau = t(:) - t(1);

    lambda = hodmd.lambda(:);
    lambdaSafe = lambda;
    lambdaSafe(abs(lambdaSafe) < eps) = eps;
    omega = log(lambdaSafe) / dt;

    [G, basisInfo] = build_modal_time_basis(tau, omega, lambda, modeIdx);

    if isempty(G)
        return;
    end

    sG = svd(G, 'econ');

    if isempty(sG) || sG(1) == 0
        condG = inf;
        rankG = 0;
    else
        tolG = max(size(G)) * eps(sG(1));
        rankG = sum(sG > tolG);

        if sG(end) <= tolG
            condG = inf;
        else
            condG = sG(1) / sG(end);
        end
    end

    C_T = pinv(G) * A_rec.';
    C = C_T.';

    A_fit = C * G.';
    fitResidual = A_rec - A_fit;

    typeList = {basisInfo.type};
    C_real = C(:, strcmp(typeList, 'real'));
    C_cos = C(:, strcmp(typeList, 'cos'));
    C_sin = C(:, strcmp(typeList, 'sin'));

    if ~isempty(C_cos) && size(C_cos, 2) == size(C_sin, 2)
        C_amp = sqrt(C_cos.^2 + C_sin.^2);
        C_phase = atan2(C_sin, C_cos);
        qCombined = sqrt(sum(C_amp.^2, 2));
    elseif ~isempty(C_real)
        C_amp = [];
        C_phase = [];
        if size(C_real, 2) == 1
            qCombined = C_real(:, 1);
        else
            qCombined = sqrt(sum(C_real.^2, 2));
        end
    else
        C_amp = [];
        C_phase = [];
        qCombined = sqrt(sum(C.^2, 2));
    end

    qSpectrum.valid = true;
    qSpectrum.G = G;
    qSpectrum.C = C;
    qSpectrum.basisInfo = basisInfo;
    qSpectrum.condG = condG;
    qSpectrum.rankG = rankG;
    qSpectrum.qCombined = qCombined;
    qSpectrum.qModes = C;
    qSpectrum.qModeReal = C;
    qSpectrum.qModeImag = zeros(size(C));
    qSpectrum.qModeAbs = abs(C);
    qSpectrum.C_real = C_real;
    qSpectrum.C_cos = C_cos;
    qSpectrum.C_sin = C_sin;
    qSpectrum.C_amp = C_amp;
    qSpectrum.C_phase = C_phase;
    qSpectrum.A_fit = A_fit;
    qSpectrum.fitResidual = fitResidual;
    qSpectrum.fitResidualNorm = norm(fitResidual, 'fro') / max(norm(A_rec, 'fro'), eps);
end

function [G, basisInfo] = build_modal_time_basis(tau, omega, lambda, modeIdx)
    tau = tau(:);
    modeIdx = modeIdx(:).';
    used = false(size(modeIdx));
    tol = 1e-6;
    c_cm_per_fs = 2.99792458e-5;

    G = [];
    basisInfo = repmat(struct( ...
        'type', '', ...
        'group', [], ...
        'sourceMode', [], ...
        'pairedMode', [], ...
        'omega', [], ...
        'growthRate_fs', [], ...
        'angularFrequency_rad_per_fs', [], ...
        'frequency_fs', [], ...
        'period_fs', [], ...
        'wavenumber_cm1', []), 0, 1);

    group = 0;

    for ii = 1:numel(modeIdx)
        if used(ii)
            continue;
        end

        idx = modeIdx(ii);
        partner = find_conjugate_partner(lambda, idx, tol);
        if isempty(partner)
            partnerPos = [];
        else
            partnerPos = find(modeIdx == partner, 1);
        end
        hasSelectedPartner = ~isempty(partner) && ~isempty(partnerPos) && partner ~= idx;

        if hasSelectedPartner
            used(ii) = true;
            used(partnerPos) = true;

            if imag(omega(idx)) >= imag(omega(partner))
                sourceMode = idx;
                pairedMode = partner;
            else
                sourceMode = partner;
                pairedMode = idx;
            end

            om = omega(sourceMode);
            group = group + 1;
            [G, basisInfo] = append_oscillatory_basis( ...
                G, basisInfo, tau, om, sourceMode, pairedMode, group, c_cm_per_fs);
            continue;
        end

        used(ii) = true;
        om = omega(idx);
        group = group + 1;

        if abs(imag(om)) > 1e-12
            [G, basisInfo] = append_oscillatory_basis( ...
                G, basisInfo, tau, om, idx, nan, group, c_cm_per_fs);
        else
            alpha = real(om);
            G(:, end+1) = exp(alpha * tau); %#ok<AGROW>
            basisInfo(end+1, 1) = make_basis_info( ...
                'real', group, idx, nan, om, c_cm_per_fs); %#ok<AGROW>
        end
    end
end

function [G, basisInfo] = append_oscillatory_basis( ...
    G, basisInfo, tau, om, sourceMode, pairedMode, group, c_cm_per_fs)
    alpha = real(om);
    Omega = abs(imag(om));

    G(:, end+1) = exp(alpha * tau) .* cos(Omega * tau);
    basisInfo(end+1, 1) = make_basis_info( ...
        'cos', group, sourceMode, pairedMode, om, c_cm_per_fs);

    G(:, end+1) = exp(alpha * tau) .* sin(Omega * tau);
    basisInfo(end+1, 1) = make_basis_info( ...
        'sin', group, sourceMode, pairedMode, om, c_cm_per_fs);
end

function info = make_basis_info(type, group, sourceMode, pairedMode, om, c_cm_per_fs)
    freq_fs = abs(imag(om)) / (2*pi);

    if freq_fs > eps
        period_fs = 1 / freq_fs;
    else
        period_fs = inf;
    end

    info = struct();
    info.type = type;
    info.group = group;
    info.sourceMode = sourceMode;
    info.pairedMode = pairedMode;
    info.omega = om;
    info.growthRate_fs = real(om);
    info.angularFrequency_rad_per_fs = imag(om);
    info.frequency_fs = freq_fs;
    info.period_fs = period_fs;
    info.wavenumber_cm1 = freq_fs / c_cm_per_fs;
end

function fI = i2_atomic_form_factor(q)
    coeff = [20.1472 4.347 18.9949 0.3814 7.513 ...
             27.766 2.2735 66.8776 4.0712];
    a = coeff([1 3 5 7]);
    b = coeff([2 4 6 8]);
    c = coeff(9);

    x = (q(:) ./ (4*pi)).^2;
    fI = sum(a .* exp(-x * b), 2) + c;
end


%% ============================================================
% PLOTTING / INSPECTION
% =============================================================



