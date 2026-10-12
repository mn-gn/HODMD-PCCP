%% step02_run_HODMD.m
% HODMD/DMD comparison for TRXSS difference signal A(q,t)
% This script compares the conventional first-delay-vector amplitudes b_j
% with globally fitted amplitudes beta_j over all delay vectors.
%
% Input .mat variables:
%   A : Nq x Nt
%   q : Nq x 1
%   t : 1 x Nt or Nt x 1, stored in seconds
%
% Pipeline:
%   1. Load A(q,t)
%   2. Crop q and time using the configured ranges/index
%   3. Convert time axis from seconds to femtoseconds for display/analysis
%   4. Optional preprocessing
%   5. Use preprocessed A directly; no SVD coordinate reduction
%   6. Delay embedding, or direct snapshots when delayOrder = 0
%   7. DMD on delayed/direct snapshots
%   8. Mode table: period, frequency, decay/growth rate
%   9. Fit both b_j and beta_j
%  10. Compare mode 2,3 pair 2D reconstructions from b_j vs beta_j

clearvars; clc; close all;
set(groot, 'DefaultFigureWindowStyle', 'docked');

%% ============================================================
% 0. USER CONFIGURATION
% =============================================================

cfg = struct();

% ---------- input ----------
projectDir = fileparts(mfilename('fullpath'));
cfg.matFile = fullfile(projectDir, 'data', 'mock', 'I2_qDeltaS_noise90_seedN1.mat');
cfg.varA = 'A';
cfg.varq = 'q';
cfg.vart = 't';

% ---------- crop ----------
cfg.cropStartIndex = 1;  % use A(:,cropStartIndex:end), t(cropStartIndex:end)
cfg.cropQRange = [];      % []: all q; [1.5 5.5]: inclusive q range (Angstrom^-1)

% ---------- time unit ----------
% Input t is in seconds, so multiply by 1e15 to display/analyze in fs.
cfg.fsPerInputTimeUnit = 1e15;

% ---------- preprocessing: q moving average only ----------
cfg.pre.smoothAlongQ = true;
cfg.pre.qSmoothWindow = 20; % number of q samples, endpoints shrink

% ---------- HODMD/DMD parameters ----------
% delay order d. Use d = 0 for standard DMD without delay embedding.
% crop 후 Nt=81이면 10~25 정도 scan.
cfg.hodmd.delayOrder = 20;

% DMD rank after delay embedding, or after direct snapshots when d = 0.
% 너무 크면 noise/growing mode가 많아짐.
cfg.hodmd.dmdRank = 9;

% ---------- b_j vs beta_j comparison ----------
% b_j uses only the first delay vector:
%   Z(:,1) ~= Phi * b
% beta_j uses all delay vectors globally:
%   Z(:,k) ~= sum_j Phi(:,j) beta_j lambda_j^(k-1)
cfg.compareBeta.enabled = true;
cfg.compareBeta.modeIdx = [2 3];
cfg.compareBeta.plotResults = true;
cfg.compareBeta.onlyCompare = false;

% ---------- restricted global beta comparison ----------
% Refit the global amplitudes after retaining only mode 1, or modes 1,2,3.
% Mode numbers refer to the ordering obtained from the all-mode beta fit.
cfg.restrictedBetaComparison.enabled = true;

% ---------- mode(2,3) coherent-amplitude diagnostic ----------
% Compare the current mode(2,3) reconstruction with the mode(2,3)
% reconstruction obtained from the matching zero-noise mock dataset.
cfg.mode23Projection.enabled = true;
cfg.mode23Projection.modeIdx = [2 3];
cfg.mode23Projection.referenceMatFile = '';
% Empty referenceMatFile automatically replaces "noise###" with "noise0"
% in cfg.matFile while preserving the seed-count suffix.

% ---------- target mode search ----------
cfg.target.period_fs = 303.24;
cfg.target.tolerance_fs = 20;

% 각 target period별로 positive-frequency mode 몇 개를 자동 선택할지
cfg.target.numPositiveModesPerTarget = 1;

% ---------- I2 mock distance comparison ----------
% These parameters match step01_generate_noise.m.
% The inverse estimate uses the target-mode oscillatory q signal and the
% linearized Debye response around R_eq.
cfg.mockI2Distance.enabled = false;
cfg.mockI2Distance.R_eq_A = 3.300;
cfg.mockI2Distance.amp_A = 0.1;
cfg.mockI2Distance.period_fs = 303.24;
cfg.mockI2Distance.tau_fs = 400;
cfg.mockI2Distance.phi_rad = pi;
cfg.mockI2Distance.R_ground_A = 2.700;
cfg.mockI2Distance.plotResults = false;

% ---------- manual mode reconstruction ----------
% 비워두면 target period 근처 mode를 자동으로 복원.
% mode 번호는 출력되는 mode table의 Mode 번호 기준.
%
% 예:
%   cfg.manualSelectedModes = [3 4];
%   cfg.manualSelectedModes = [7 8 9 10];
cfg.manualSelectedModes = [];

% manualSelectedModes가 있을 때 conjugate pair도 자동 포함할지
cfg.includeConjugatePairsForManual = true;

% ---------- plotting ----------
cfg.plot.maxModesToShow = 40;
cfg.plot.useSymmetricColorLimit = true;
cfg.plot.percentileForCLim = 99.5;
cfg.plot.maxFreqToPlot = 0.02; % cycles/fs

% ---------- save ----------
% Keep this false so only the global-beta structure-fit MAT file is saved.
cfg.saveResults = false;
cfg.outputDir = fullfile(fileparts(mfilename('fullpath')), 'results');
% All beta_j outputs are collected in the method-local results folder.
cfg.outputFile = [];
% If enabled manually, outputFile uses <inputMatName>_hodmdResults.mat.

% ---------- q-spectrum save ----------
% true로 바꾸면 manual/automatic reconstruction에서 plot되는
% combined real q-spectrum과 q를 .mat 파일로 저장함.
cfg.qSpectrum.saveMat = false;
cfg.qSpectrum.outputDir = [];
% Empty outputDir uses cfg.outputDir.
cfg.qSpectrum.filePrefix = [];
% Empty filePrefix uses <inputMatName>_globalBeta_qSpec.

% true이면 qModes도 같이 저장함.
% qModes(:,j)는 선택된 각 mode의 complex individual q-spectrum.
cfg.qSpectrum.saveIndividualModes = false;

% Sine transform of signed modal q coefficients (q in Angstrom^-1).
cfg.rSpectrum.enabled = true;
cfg.rSpectrum.inputIsQWeighted = true; % true: A = q*DeltaS; false: A = DeltaS
cfg.rSpectrum.range_A = [0 10];
cfg.rSpectrum.step_A = 0.01;
cfg.rSpectrum.evaluationTime_fs = []; % []: analysis start t(1); scalar: absolute fs
% F(r) = (2/pi) integral C(q)*sin(q*r) dq for q-weighted input.
% Uses the measured q range without extrapolation or an extra window.

% ---------- selected 2D mode-pair reconstruction save ----------
% Saves selected mode reconstruction A(q,t), q, t, period, and damping time.
% Empty outputDir uses cfg.outputDir.
cfg.modePairSave.enabled = false;
cfg.modePairSave.modeIdx = [2 3];
cfg.modePairSave.outputDir = [];
cfg.modePairSave.outputFile = [];
% Empty outputFile uses <inputMatName>_mode2_3.mat.

% ---------- structure-fitting 2D component save ----------
% Saves de-Hankelized/reconstructed 2D components A(q,t) for structure fitting.
% This avoids interpreting a single extracted q-spectrum as the physical signal.
% Empty outputDir uses cfg.outputDir.
cfg.structureFitSave.enabled = false;
cfg.structureFitSave.trendModeIdx = 1;
cfg.structureFitSave.oscillationModeIdx = [4 5];
cfg.structureFitSave.outputDir = [];
cfg.structureFitSave.outputFile = [];
% Empty outputFile uses <inputMatName>_globalBeta_structureFit.mat.


%% ============================================================
% 1. LOAD DATA
% =============================================================

S = load(cfg.matFile);

cfg.inputMetadata = struct('hasFWHM', false, 'FWHMUnit', 's');
if isfield(S, 'FWHM')
    FWHM = S.FWHM;
    assert(isnumeric(FWHM) && isreal(FWHM) && ~isempty(FWHM) && ...
        all(isfinite(FWHM(:))) && all(FWHM(:) > 0), ...
        'Input variable FWHM must contain finite positive numeric value(s) in seconds.');

    cfg.inputMetadata.hasFWHM = true;
    cfg.inputMetadata.FWHM = FWHM;

    FWHM_fs_display = double(FWHM(:)) .* 1e15;
    if isscalar(FWHM)
        fprintf('Loaded input FWHM metadata: %.6g s = %.6g fs\n', ...
            double(FWHM), FWHM_fs_display);
    else
        fprintf(['Loaded input FWHM metadata: %d values in seconds ', ...
            '(%.6g to %.6g fs)\n'], ...
            numel(FWHM), min(FWHM_fs_display), max(FWHM_fs_display));
    end
end


q = S.(cfg.varq);
A = S.(cfg.varA);
t = S.(cfg.vart);

q = q(:);
t = t(:).';

if ndims(A) > 3
    error('A must be 2D or 3D. For 3D input, use A as Nq x Nt x Nseries.');
end

if size(A,1) == numel(q) && size(A,2) == numel(t)
    % A is already Nq x Nt, or Nq x Nt x Nseries.
elseif size(A,1) == numel(t) && size(A,2) == numel(q)
    if ismatrix(A)
        A = A.';
    else
        A = permute(A, [2 1 3]);
    end
    fprintf('Transposed A so that rows match q and columns match t.\n');
else
    error('Cannot orient A. size(A) must match q x t or t x q, with optional 3rd axis.');
end

assert(all(isfinite(A(:))), 'A contains NaN or Inf values.');

if any(~isfinite(q)) || any(~isfinite(t))
    error('q and t must contain only finite values.');
end

q_input = q;
[A, q, qCropInfo] = crop_q_range(A, q, cfg.cropQRange);
if qCropInfo.applied
    fprintf('Cropped q to [%.6g, %.6g] Angstrom^-1 (%d of %d points retained).\n', ...
        min(q), max(q), numel(q), numel(q_input));
end
if isempty(cfg.fsPerInputTimeUnit) || ~isscalar(cfg.fsPerInputTimeUnit) || ...
        ~isfinite(cfg.fsPerInputTimeUnit) || cfg.fsPerInputTimeUnit <= 0
    error('cfg.fsPerInputTimeUnit must be a finite positive scalar.');
end

% sort time axis
[t, sortIdx] = sort(t, 'ascend');
if ndims(A) == 3
    A = A(:, sortIdx, :);
else
    A = A(:, sortIdx);
end

% crop
idx0 = cfg.cropStartIndex;

assert(idx0 >= 1 && idx0 <= numel(t), ...
    'cfg.cropStartIndex must be between 1 and Nt.');

t_seconds = t(idx0:end);
t = t_seconds * cfg.fsPerInputTimeUnit;
if ndims(A) == 3
    A = A(:, idx0:end, :);
else
    A = A(:, idx0:end);
end

if ndims(A) == 3
    [Nq, Nt, Nseries] = size(A);

    assert(numel(q) == Nq, 'q length does not match size(A,1).');
    assert(numel(t) == Nt, 't length does not match size(A,2).');

    dtVec = diff(t);
    dt = median(dtVec);

    if max(abs(dtVec - dt)) > 1e-6 * max(abs(dt), eps)
        warning('t is not perfectly equally spaced. HODMD assumes approximately uniform time spacing.');
    end

    cfg.hodmd.delayOrder = round(cfg.hodmd.delayOrder);

    assert(cfg.hodmd.delayOrder == 0 || ...
        (cfg.hodmd.delayOrder >= 2 && cfg.hodmd.delayOrder <= Nt-2), ...
        ['cfg.hodmd.delayOrder must be 0 for standard DMD, ' ...
         'or satisfy 2 <= delayOrder <= Nt-2 for HODMD.']);

    fprintf('\nLoaded and cropped 3D data:\n');
    fprintf('  A size after crop = %d x %d x %d\n', Nq, Nt, Nseries);
    fprintf('  t range           = %.6g to %.6g s\n', t_seconds(1), t_seconds(end));
    fprintf('  displayed t range = %.6g to %.6g fs\n', t(1), t(end));
    fprintf('  dt                = %.6g fs\n', dt);
    fprintf('\n');

    batchResults = run_hodmd_batch_3d(A, q, t, t_seconds, dt, cfg);
    batchResults.qCropInfo = qCropInfo;

    if cfg.saveResults
        resultsOutputFile = output_file_in_dir( ...
            cfg, cfg.outputFile, 'hodmdResults', global_beta_output_dir(cfg));

        save(resultsOutputFile, ...
            'A', 'q', 'q_input', 'qCropInfo', ...
            't', 't_seconds', 'cfg', 'batchResults', ...
            '-v7.3');
        append_input_fwhm_if_available(resultsOutputFile, cfg);

        fprintf('\nSaved 3D batch results to %s\n', resultsOutputFile);
    end

    return;
end

[Nq, Nt] = size(A);

assert(numel(q) == Nq, 'q length does not match size(A,1).');
assert(numel(t) == Nt, 't length does not match size(A,2).');

dtVec = diff(t);
dt = median(dtVec);

if max(abs(dtVec - dt)) > 1e-6 * max(abs(dt), eps)
    warning('t is not perfectly equally spaced. HODMD assumes approximately uniform time spacing.');
end

fprintf('\nLoaded and cropped data:\n');
fprintf('  A size after crop = %d x %d\n', Nq, Nt);
fprintf('  t range           = %.6g to %.6g s\n', t_seconds(1), t_seconds(end));
fprintf('  displayed t range = %.6g to %.6g fs\n', t(1), t(end));
fprintf('  dt                = %.6g fs\n', dt);
fprintf('\n');

cfg.hodmd.delayOrder = round(cfg.hodmd.delayOrder);

assert(cfg.hodmd.delayOrder == 0 || ...
    (cfg.hodmd.delayOrder >= 2 && cfg.hodmd.delayOrder <= Nt-2), ...
    ['cfg.hodmd.delayOrder must be 0 for standard DMD, ' ...
     'or satisfy 2 <= delayOrder <= Nt-2 for HODMD.']);


%% ============================================================
% 2. PREPROCESSING
% =============================================================

A_original = A;
A_work = preprocess_hodmd_matrix(A, q, cfg);
A_scale = ones(Nq, 1); % retained in reconstruction interfaces; no normalization

plot_qt_map(t, q, A_original, 'Original cropped A(q,t)', cfg);
plot_qt_map(t, q, A_work, 'Input to HODMD after preprocessing', cfg);

% HODMD acts directly on the smoothed q-space data.
Uq = eye(Nq);
B = A_work;
rState = Nq;

% Delay embedding; the DMD-input Z1 singular spectrum is retained.
d = cfg.hodmd.delayOrder;

if d == 0
    Z = B;
else
    Z = make_delay_matrix(B, d);
end

[nAug, K] = size(Z);

assert(K >= 2, 'DMD/HODMD requires at least two snapshots after preprocessing.');

if d == 0
    fprintf('Standard DMD input without delay embedding:\n');
    fprintf('  delay order d = 0\n');
    fprintf('  state dimension = %d\n', rState);
    fprintf('  number of snapshots K = %d\n\n', K);
else
    fprintf('Delay embedding:\n');
    fprintf('  delay order d = %d\n', d);
    fprintf('  state dimension = %d\n', rState);
    fprintf('  augmented dimension = %d\n', nAug);
    fprintf('  number of delay snapshots K = %d\n\n', K);
end

Z1 = Z(:, 1:end-1);
Z2 = Z(:, 2:end);


%% ============================================================
% 5. DMD ON DELAYED OR DIRECT SNAPSHOTS
% =============================================================

[hodmd, dmdInfo] = exact_dmd(Z1, Z2, cfg.hodmd.dmdRank, dt);
dmdInfo.delayOrder = d;
dmdInfo.isDelayEmbedded = d > 0;

if d == 0
    fprintf('DMD done:\n');
else
    fprintf('HODMD done:\n');
end
fprintf('  requested DMD rank = %d\n', cfg.hodmd.dmdRank);
fprintf('  actual DMD rank    = %d\n\n', numel(hodmd.lambda));

% Compute modal amplitudes by both conventions.
[betaGlobal, bFirstSnapshot, betaFitInfo] = ...
    fit_global_dmd_amplitudes(hodmd.Phi, hodmd.lambda, Z);

hodmd_b = hodmd;
hodmd_b.b = bFirstSnapshot;
hodmd_b.b_firstSnapshot = bFirstSnapshot;
hodmd_b.beta = betaGlobal;
hodmd_b.betaFitInfo = betaFitInfo;
hodmd_b.amplitudeMethod = 'first_delay_vector_b';

hodmd_beta = hodmd;
hodmd_beta.b = betaGlobal;
hodmd_beta.beta = betaGlobal;
hodmd_beta.b_firstSnapshot = bFirstSnapshot;
hodmd_beta.betaFitInfo = betaFitInfo;
hodmd_beta.amplitudeMethod = 'global_delay_vector_beta';

% Continue the standard mode table/diagnostics with the global beta
% amplitudes, but keep the b_j version for direct comparison.
hodmd = hodmd_beta;

fprintf('  first-delay-vector amplitude fit rel. error over Z = %.6g\n', ...
    betaFitInfo.relErr_firstSnapshotB);
fprintf('  global beta amplitude fit rel. error over Z       = %.6g\n', ...
    betaFitInfo.relErr_globalBeta);
fprintf('  global beta LS rank / condition                  = %d / %.6g\n\n', ...
    betaFitInfo.rankX, betaFitInfo.condX);

% Sort modes by global-beta modal amplitude and apply the same ordering to
% the b_j copy. This keeps mode numbers aligned between the two methods.
[hodmd, betaSortIdx] = sort_hodmd_modes_by_amplitude(hodmd);
hodmd_beta = hodmd;
hodmd_b = apply_hodmd_sort_index(hodmd_b, betaSortIdx);

restrictedBetaComparison = struct();
if isfield(cfg, 'restrictedBetaComparison') && ...
        cfg.restrictedBetaComparison.enabled
    restrictedBetaComparison = compare_restricted_global_beta(hodmd, Z);
    print_restricted_global_beta_comparison(restrictedBetaComparison);
end

modeTable = build_hodmd_mode_table(hodmd, dt);

if d == 0
    fprintf('DMD mode table, sorted by modal amplitude:\n');
else
    fprintf('HODMD mode table, sorted by modal amplitude:\n');
end
disp(modeTable(1:min(cfg.plot.maxModesToShow, height(modeTable)), :));

mode23Projection = struct();
if isfield(cfg, 'mode23Projection') && cfg.mode23Projection.enabled
    mode23Reference = prepare_mode23_projection_reference(q, t, dt, cfg);

    if mode23Reference.valid
        projectionModeIdx = unique_valid_indices( ...
            cfg.mode23Projection.modeIdx, numel(hodmd.lambda));
        B_mode23Projection = reconstruct_hodmd_modes( ...
            hodmd, projectionModeIdx, rState, d, Nt);
        A_mode23Projection = Uq * B_mode23Projection;
        A_mode23Projection = undo_component_scale( ...
            A_mode23Projection, A_scale, cfg);

        mode23Projection = compute_mode23_projection_diagnostics( ...
            A_mode23Projection, mode23Reference.A_mode23, ...
            projectionModeIdx, mode23Reference);
        print_mode23_projection_diagnostics(mode23Projection);
    end
end

plot_hodmd_diagnostics(hodmd, dmdInfo, modeTable, cfg);

betaComparison = struct();
if isfield(cfg, 'compareBeta') && cfg.compareBeta.enabled
    betaComparison = compare_b_beta_mode_pair( ...
        hodmd_b, hodmd_beta, q, t, Uq, A_scale, cfg, rState, d, Nt);

    if isfield(cfg.compareBeta, 'onlyCompare') && cfg.compareBeta.onlyCompare
        fprintf('\nComparison-only mode is enabled. No result files were saved.\n');
        return;
    end
end


%% ============================================================
% 6. TARGET PERIOD MODE SEARCH
% =============================================================

targetGroups = find_target_hodmd_modes(hodmd, modeTable, cfg.target);

fprintf('\nTarget period candidate modes:\n');
c_cm_per_s = 2.99792458e10;  % speed of light [cm/s]
for i = 1:numel(cfg.target.period_fs)
    target_cm1 = 1e15 / (cfg.target.period_fs(i) * c_cm_per_s);
    fprintf('  target %.1f fs (%.2f cm^-1) -> modes %s\n', ...
        cfg.target.period_fs(i), target_cm1, mat2str(targetGroups{i}));
end


%% ============================================================
% 7. RECONSTRUCT MODES
% =============================================================

recon = struct();

modeComponentRecon = reconstruct_all_mode_components( ...
    hodmd, Uq, A_scale, cfg, rState, d, Nt, modeTable);
A_modeComponents = modeComponentRecon.A;
modeComponentTable = modeComponentRecon.table;

recon.modeComponents = modeComponentRecon;

fprintf('\nDynamic modes ranked by summed eigencomponent amplitude:\n');
disp(modeComponentTable);

plot_all_mode_component_reconstructions(t, q, modeComponentRecon, cfg);

structureFitSave = save_structure_fit_components_if_requested( ...
    q, t, hodmd, Uq, A_scale, cfg, rState, d, Nt, modeTable);
recon.structureFitSave = structureFitSave;

selectedModePairSave = save_selected_mode_pair_reconstruction_if_requested( ...
    q, t, hodmd, Uq, A_scale, cfg, rState, d, Nt, modeTable);
recon.selectedModePairSave = selectedModePairSave;

targetModeIdx = [];
for i = 1:numel(targetGroups)
    targetModeIdx = [targetModeIdx, targetGroups{i}]; %#ok<AGROW>
end
targetModeIdx = unique_valid_indices(targetModeIdx, numel(hodmd.lambda));

B_targetComparison = reconstruct_hodmd_modes(hodmd, targetModeIdx, rState, d, Nt);
A_targetComparison = Uq * B_targetComparison;
A_targetComparison = undo_component_scale(A_targetComparison, A_scale, cfg);

recon.targetComparison.targetModeIdx = targetModeIdx;
recon.targetComparison.A = A_targetComparison;

if isfield(cfg, 'mockI2Distance') && cfg.mockI2Distance.enabled
    if isempty(targetModeIdx)
        warning('I2 distance comparison skipped because no target modes were selected.');
    else
        i2DistanceComparison = estimate_i2_distance_from_target_modes( ...
            t, q, A_targetComparison, targetModeIdx, cfg);
        recon.targetComparison.i2Distance = i2DistanceComparison;

        if isfield(cfg.mockI2Distance, 'plotResults') && ...
                cfg.mockI2Distance.plotResults
            plot_i2_distance_comparison(i2DistanceComparison);
        end
    end
end

% Manual reconstruction
rSpectrumFields = {};
if ~isempty(cfg.manualSelectedModes)
    modeIdx = cfg.manualSelectedModes(:).';

    if cfg.includeConjugatePairsForManual
        modeIdx = include_conjugate_pairs(hodmd.lambda, modeIdx);
    end

    modeIdx = unique_valid_indices(modeIdx, numel(hodmd.lambda));

    fprintf('\nManual selected modes for reconstruction:\n');
    disp(modeIdx);

    B_manual = reconstruct_hodmd_modes(hodmd, modeIdx, rState, d, Nt);
    A_manual = Uq * B_manual;
    A_manual = undo_component_scale(A_manual, A_scale, cfg);

    recon.manual.modeIdx = modeIdx;
    recon.manual.A = A_manual;

    ttl = sprintf('HODMD manual selected modes %s', mat2str(modeIdx));
    plot_qt_map(t, q, A_manual, ttl, cfg);
    inspect_reconstructed_signal(A_manual, t, cfg, ttl);
    recon.manual.qSpectrum = plot_hodmd_q_spectrum( ...
        q, t, Uq, hodmd, modeIdx, rState, d, A_scale, cfg, ttl);
    rSpectrumFields{end+1} = 'manual';
end

% Automatic target reconstructions
for i = 1:numel(cfg.target.period_fs)
    target_fs = cfg.target.period_fs(i);
    modeIdx = targetGroups{i};

    if isempty(modeIdx)
        continue;
    end

    B_target = reconstruct_hodmd_modes(hodmd, modeIdx, rState, d, Nt);
    A_target = Uq * B_target;
    A_target = undo_component_scale(A_target, A_scale, cfg);

    fieldName = sprintf('target_%dfs', round(target_fs));
    recon.(fieldName).modeIdx = modeIdx;
    recon.(fieldName).A = A_target;

    ttl = sprintf('HODMD reconstructed target %.0f fs, modes %s', ...
        target_fs, mat2str(modeIdx));

    plot_qt_map(t, q, A_target, ttl, cfg);
    inspect_reconstructed_signal(A_target, t, cfg, ttl);
    recon.(fieldName).qSpectrum = plot_hodmd_q_spectrum( ...
        q, t, Uq, hodmd, modeIdx, rState, d, A_scale, cfg, ttl);
    rSpectrumFields{end+1} = fieldName; %#ok<SAGROW>
end

% Final figures: transform each signed temporal-basis coefficient separately.
if cfg.rSpectrum.enabled
    for i = 1:numel(rSpectrumFields)
        fieldName = rSpectrumFields{i};
        recon.(fieldName).rSpectrum = plot_hodmd_r_spectrum( ...
            recon.(fieldName).qSpectrum, cfg.rSpectrum, fieldName);
    end
end


%% ============================================================
% 8. OPTIONAL SAVE
% =============================================================

if cfg.saveResults
    resultsOutputFile = output_file_in_dir( ...
        cfg, cfg.outputFile, 'hodmdResults', global_beta_output_dir(cfg));

    save(resultsOutputFile, ...
        'A_original', 'A_work', ...
        'q', 'q_input', 'qCropInfo', 't', 't_seconds', 'cfg', ...
        'Uq', 'B', ...
        'Z', 'hodmd', 'hodmd_b', 'hodmd_beta', 'betaComparison', ...
        'restrictedBetaComparison', ...
        'mode23Projection', ...
        'dmdInfo', 'modeTable', 'targetGroups', ...
        'modeComponentRecon', 'A_modeComponents', 'modeComponentTable', ...
        'structureFitSave', 'selectedModePairSave', 'recon', ...
        '-v7.3');
    append_input_fwhm_if_available(resultsOutputFile, cfg);

    fprintf('\nSaved results to %s\n', resultsOutputFile);
end


%% ============================================================
% LOCAL FUNCTIONS
% =============================================================

function batchResults = run_hodmd_batch_3d(A3, q, t, ~, dt, cfg)
    [Nq, Nt, nSeries] = size(A3);
    % Select the complete oscillatory pair nearest 110 cm^-1, without a gate.
    cfg.batchTargetWavenumber_cm1 = 110;
    targetPeriods = 1 / (2.99792458e-5 * cfg.batchTargetWavenumber_cm1);
    nTargets = 1;
    saveModePair = isfield(cfg, 'modePairSave') && ...
        isfield(cfg.modePairSave, 'enabled') && cfg.modePairSave.enabled;
    saveStructureFit = is_structure_fit_save_enabled(cfg);
    runMode23Projection = isfield(cfg, 'mode23Projection') && ...
        cfg.mode23Projection.enabled;

    mode23Reference = struct('valid', false);
    if runMode23Projection
        mode23Reference = prepare_mode23_projection_reference(q, t, dt, cfg);
    end

    Frequency_fs = nan(nSeries, nTargets);
    Wavenumber_cm1 = nan(nSeries, nTargets);
    Period_fs = nan(nSeries, nTargets);
    GrowthRate_fs = nan(nSeries, nTargets);
    DecayTime_fs = nan(nSeries, nTargets);
    Amplitude = nan(nSeries, nTargets);
    ModeIndex = nan(nSeries, nTargets);
    TargetModeGroup = cell(nSeries, nTargets);

    success = false(nSeries, 1);
    errorMessage = cell(nSeries, 1);
    perRunModeTable = cell(nSeries, 1);
    perRunDynamicModeTable = cell(nSeries, 1);
    perRunTargetGroups = cell(nSeries, 1);
    perRunBetaFitInfo = cell(nSeries, 1);
    perRunRestrictedBetaComparison = cell(nSeries, 1);

    mode23ProjectionPerRun = cell(nSeries, 1);
    mode23TargetNorm = nan(nSeries, 1);
    mode23ParallelNorm = nan(nSeries, 1);
    mode23OrthogonalNorm = nan(nSeries, 1);
    mode23ProjectionCoefficient = nan(nSeries, 1);
    mode23ParallelEnergyFraction = nan(nSeries, 1);
    mode23OrthogonalEnergyFraction = nan(nSeries, 1);
    mode23CosineSimilarity = nan(nSeries, 1);

    if saveModePair
        nSaveModes = numel(cfg.modePairSave.modeIdx);
        nSaveBasisMax = max(1, 2*nSaveModes);
        A_mode23_all = nan(size(A3));
        selectedModeIdx = nan(nSeries, numel(cfg.modePairSave.modeIdx));
        representativeMode = nan(nSeries, 1);
        selectedPeriod_fs = nan(nSeries, 1);
        selectedDecayTime_fs = nan(nSeries, 1);
        selectedWavenumber_cm1 = nan(nSeries, 1);
        selectedFrequency_fs = nan(nSeries, 1);
        selectedGrowthRate_fs = nan(nSeries, 1);
        qSpectrumCombined_all = nan(Nq, nSeries);
        qSpectrumModes_all = nan(Nq, nSaveBasisMax, nSeries);
        qSpectrumBasisCoefficients_all = nan(Nq, nSaveBasisMax, nSeries);
        qSpectrumModeAbs_all = nan(Nq, nSaveBasisMax, nSeries);
        qSpectrumModeReal_all = nan(Nq, nSaveBasisMax, nSeries);
        qSpectrumModeImag_all = nan(Nq, nSaveBasisMax, nSeries);
        qSpectrumCos_all = nan(Nq, nSaveModes, nSeries);
        qSpectrumSin_all = nan(Nq, nSaveModes, nSeries);
        qSpectrumAmp_all = nan(Nq, nSaveModes, nSeries);
        qSpectrumPhase_all = nan(Nq, nSaveModes, nSeries);
        qSpectrumNumBasis = nan(nSeries, 1);
        qSpectrumCondG = nan(nSeries, 1);
        qSpectrumRankG = nan(nSeries, 1);
        qSpectrumBasisInfo = cell(nSeries, 1);
    end

    if saveStructureFit
        nTrendModes = numel(cfg.structureFitSave.trendModeIdx);
        nOscillationModes = numel(cfg.structureFitSave.oscillationModeIdx);

        A_trend_all = nan(Nq, Nt, nSeries);
        A_oscillation_all = nan(Nq, Nt, nSeries);
        A_trend_plus_oscillation_all = nan(Nq, Nt, nSeries);
        trendModeIdx_all = nan(nSeries, nTrendModes);
        oscillationModeIdx_all = nan(nSeries, nOscillationModes);
        trendRepresentativeMode = nan(nSeries, 1);
        oscillationRepresentativeMode = nan(nSeries, 1);
        trendPeriod_fs = nan(nSeries, 1);
        trendDecayTime_fs = nan(nSeries, 1);
        oscillationPeriod_fs = nan(nSeries, 1);
        oscillationDecayTime_fs = nan(nSeries, 1);
        oscillationWavenumber_cm1 = nan(nSeries, 1);
    end

    fprintf('3D input detected. Running batch HODMD over %d slices.\n', nSeries);

    for iSeries = 1:nSeries
        if iSeries == 1 || iSeries == nSeries || mod(iSeries, 10) == 0
            fprintf('  analyzing slice %d / %d\n', iSeries, nSeries);
        end

        try
            singleResult = analyze_single_hodmd_target_metrics(A3(:, :, iSeries), q, t, dt, cfg);
            success(iSeries) = true;
            perRunModeTable{iSeries} = singleResult.modeTable;
            perRunDynamicModeTable{iSeries} = singleResult.dynamicModeTable;
            perRunTargetGroups{iSeries} = singleResult.targetGroups;
            perRunBetaFitInfo{iSeries} = singleResult.betaFitInfo;
            perRunRestrictedBetaComparison{iSeries} = ...
                singleResult.restrictedBetaComparison;

            if runMode23Projection && mode23Reference.valid && ...
                    ~isempty(singleResult.projectionModePairA)
                projectionResult = compute_mode23_projection_diagnostics( ...
                    singleResult.projectionModePairA, ...
                    mode23Reference.A_mode23, ...
                    singleResult.projectionModePairIndices, ...
                    mode23Reference);
                mode23ProjectionPerRun{iSeries} = projectionResult;
                mode23TargetNorm(iSeries) = projectionResult.targetNorm;
                mode23ParallelNorm(iSeries) = projectionResult.parallelNorm;
                mode23OrthogonalNorm(iSeries) = projectionResult.orthogonalNorm;
                mode23ProjectionCoefficient(iSeries) = ...
                    projectionResult.projectionCoefficient;
                mode23ParallelEnergyFraction(iSeries) = ...
                    projectionResult.parallelEnergyFraction;
                mode23OrthogonalEnergyFraction(iSeries) = ...
                    projectionResult.orthogonalEnergyFraction;
                mode23CosineSimilarity(iSeries) = ...
                    projectionResult.cosineSimilarity;
            end

            if saveStructureFit && singleResult.structureFitComponents.valid
                sf = singleResult.structureFitComponents;
                nTrend = numel(sf.trendModeIdx);
                nOsc = numel(sf.oscillationModeIdx);

                A_trend_all(:, :, iSeries) = sf.A_trend;
                A_oscillation_all(:, :, iSeries) = sf.A_oscillation;
                A_trend_plus_oscillation_all(:, :, iSeries) = sf.A_trend_plus_oscillation;
                trendModeIdx_all(iSeries, 1:nTrend) = sf.trendModeIdx;
                oscillationModeIdx_all(iSeries, 1:nOsc) = sf.oscillationModeIdx;
                trendRepresentativeMode(iSeries) = sf.trendRepresentativeMode;
                oscillationRepresentativeMode(iSeries) = sf.oscillationRepresentativeMode;
                trendPeriod_fs(iSeries) = sf.trendPeriod_fs;
                trendDecayTime_fs(iSeries) = sf.trendDecayTime_fs;
                oscillationPeriod_fs(iSeries) = sf.oscillationPeriod_fs;
                oscillationDecayTime_fs(iSeries) = sf.oscillationDecayTime_fs;
                oscillationWavenumber_cm1(iSeries) = sf.oscillationWavenumber_cm1;
            end

            if saveModePair && singleResult.selectedModePair.valid
                pairResult = singleResult.selectedModePair;
                nIdx = numel(pairResult.modeIdx);

                A_mode23_all(:, :, iSeries) = pairResult.A;
                selectedModeIdx(iSeries, 1:nIdx) = pairResult.modeIdx;
                representativeMode(iSeries) = pairResult.representativeMode;
                selectedPeriod_fs(iSeries) = pairResult.period_fs;
                selectedDecayTime_fs(iSeries) = pairResult.decayTime_fs;
                selectedWavenumber_cm1(iSeries) = pairResult.wavenumber_cm1;
                selectedFrequency_fs(iSeries) = pairResult.frequency_fs;
                selectedGrowthRate_fs(iSeries) = pairResult.growthRate_fs;

                if isfield(pairResult, 'qSpectrum') && pairResult.qSpectrum.valid
                    qSpec = pairResult.qSpectrum;
                    nBasis = size(qSpec.C, 2);
                    nBasisStore = min(nBasis, size(qSpectrumModes_all, 2));

                    qSpectrumCombined_all(:, iSeries) = qSpec.qCombined;
                    qSpectrumModes_all(:, 1:nBasisStore, iSeries) = qSpec.qModes(:, 1:nBasisStore);
                    qSpectrumBasisCoefficients_all(:, 1:nBasisStore, iSeries) = qSpec.C(:, 1:nBasisStore);
                    qSpectrumModeAbs_all(:, 1:nBasisStore, iSeries) = qSpec.qModeAbs(:, 1:nBasisStore);
                    qSpectrumModeReal_all(:, 1:nBasisStore, iSeries) = qSpec.qModeReal(:, 1:nBasisStore);
                    qSpectrumModeImag_all(:, 1:nBasisStore, iSeries) = qSpec.qModeImag(:, 1:nBasisStore);
                    qSpectrumNumBasis(iSeries) = nBasis;
                    qSpectrumCondG(iSeries) = qSpec.condG;
                    qSpectrumRankG(iSeries) = qSpec.rankG;
                    qSpectrumBasisInfo{iSeries} = qSpec.basisInfo;

                    if ~isempty(qSpec.C_cos)
                        nPairStore = min(size(qSpec.C_cos, 2), size(qSpectrumCos_all, 2));
                        qSpectrumCos_all(:, 1:nPairStore, iSeries) = qSpec.C_cos(:, 1:nPairStore);
                        qSpectrumSin_all(:, 1:nPairStore, iSeries) = qSpec.C_sin(:, 1:nPairStore);
                        qSpectrumAmp_all(:, 1:nPairStore, iSeries) = qSpec.C_amp(:, 1:nPairStore);
                        qSpectrumPhase_all(:, 1:nPairStore, iSeries) = qSpec.C_phase(:, 1:nPairStore);
                    end
                end
            end

            for iTarget = 1:nTargets
                metric = singleResult.targetMetrics(iTarget);

                Frequency_fs(iSeries, iTarget) = metric.frequency_fs;
                Wavenumber_cm1(iSeries, iTarget) = metric.wavenumber_cm1;
                Period_fs(iSeries, iTarget) = metric.period_fs;
                GrowthRate_fs(iSeries, iTarget) = metric.growthRate_fs;
                DecayTime_fs(iSeries, iTarget) = metric.decayTime_fs;
                Amplitude(iSeries, iTarget) = metric.amplitude;
                ModeIndex(iSeries, iTarget) = metric.modeIndex;
                TargetModeGroup{iSeries, iTarget} = metric.modeGroupText;
            end
        catch ME
            errorMessage{iSeries} = ME.message;
            warning('Batch HODMD failed for slice %d: %s', iSeries, ME.message);
        end
    end

    summaryTable = build_batch_target_summary_table( ...
        targetPeriods, Frequency_fs, Wavenumber_cm1, Period_fs, ...
        GrowthRate_fs, DecayTime_fs, Amplitude);
    summaryTable.TargetPeriod_fs = [];
    summaryTable.TargetWavenumber_cm1 = cfg.batchTargetWavenumber_cm1;
    summaryTable.NSelected = sum(~cellfun(@isempty, perRunTargetGroups) & ...
        cellfun(@(g) ~isempty(g) && ~isempty(g{1}), perRunTargetGroups));

    batchResults = struct();
    batchResults.targetPeriods_fs = targetPeriods;
    batchResults.selectionMethod = 'nearest_wavenumber_complex_pair_without_tolerance';
    batchResults.targetWavenumber_cm1 = cfg.batchTargetWavenumber_cm1;
    batchResults.amplitudeDefinition = 'sum_j(abs(beta_j)*norm(Phi_j,2))';
    batchResults.perRunDynamicModeTable = perRunDynamicModeTable;
    batchResults.nSeries = nSeries;
    batchResults.success = success;
    batchResults.errorMessage = errorMessage;
    batchResults.Frequency_fs = Frequency_fs;
    batchResults.Wavenumber_cm1 = Wavenumber_cm1;
    batchResults.Period_fs = Period_fs;
    batchResults.GrowthRate_fs = GrowthRate_fs;
    batchResults.DecayTime_fs = DecayTime_fs;
    batchResults.Amplitude = Amplitude;
    batchResults.ModeIndex = ModeIndex;
    batchResults.TargetModeGroup = TargetModeGroup;
    batchResults.perRunModeTable = perRunModeTable;
    batchResults.perRunTargetGroups = perRunTargetGroups;
    batchResults.perRunBetaFitInfo = perRunBetaFitInfo;
    batchResults.perRunRestrictedBetaComparison = ...
        perRunRestrictedBetaComparison;
    batchResults.summaryTable = summaryTable;
    batchResults.truth = batch_truth_values(cfg);

    if runMode23Projection && mode23Reference.valid
        mode23ProjectionBatch = struct();
        mode23ProjectionBatch.reference = mode23Reference;
        mode23ProjectionBatch.perRun = mode23ProjectionPerRun;
        mode23ProjectionBatch.targetNorm = mode23TargetNorm;
        mode23ProjectionBatch.parallelNorm = mode23ParallelNorm;
        mode23ProjectionBatch.orthogonalNorm = mode23OrthogonalNorm;
        mode23ProjectionBatch.projectionCoefficient = ...
            mode23ProjectionCoefficient;
        mode23ProjectionBatch.parallelEnergyFraction = ...
            mode23ParallelEnergyFraction;
        mode23ProjectionBatch.orthogonalEnergyFraction = ...
            mode23OrthogonalEnergyFraction;
        mode23ProjectionBatch.cosineSimilarity = mode23CosineSimilarity;
        mode23ProjectionBatch.summaryTable = ...
            build_mode23_projection_batch_summary(mode23ProjectionBatch);
        batchResults.mode23Projection = mode23ProjectionBatch;
    else
        batchResults.mode23Projection = struct();
    end

    if saveModePair
        selectedModePairBatch = struct();
        selectedModePairBatch.A_mode23_all = A_mode23_all;
        selectedModePairBatch.A_selectedModes_all = A_mode23_all;
        selectedModePairBatch.q = q;
        selectedModePairBatch.t = t;
        selectedModePairBatch.selectedModeIdx = selectedModeIdx;
        selectedModePairBatch.representativeMode = representativeMode;
        selectedModePairBatch.period_fs = selectedPeriod_fs;
        selectedModePairBatch.decayTime_fs = selectedDecayTime_fs;
        selectedModePairBatch.wavenumber_cm1 = selectedWavenumber_cm1;
        selectedModePairBatch.frequency_fs = selectedFrequency_fs;
        selectedModePairBatch.growthRate_fs = selectedGrowthRate_fs;
        selectedModePairBatch.qSpectrumCombined_all = qSpectrumCombined_all;
        selectedModePairBatch.qSpectrumModes_all = qSpectrumModes_all;
        selectedModePairBatch.qSpectrumBasisCoefficients_all = qSpectrumBasisCoefficients_all;
        selectedModePairBatch.qSpectrumModeAbs_all = qSpectrumModeAbs_all;
        selectedModePairBatch.qSpectrumModeReal_all = qSpectrumModeReal_all;
        selectedModePairBatch.qSpectrumModeImag_all = qSpectrumModeImag_all;
        selectedModePairBatch.qSpectrumCos_all = qSpectrumCos_all;
        selectedModePairBatch.qSpectrumSin_all = qSpectrumSin_all;
        selectedModePairBatch.qSpectrumAmp_all = qSpectrumAmp_all;
        selectedModePairBatch.qSpectrumPhase_all = qSpectrumPhase_all;
        selectedModePairBatch.qSpectrumNumBasis = qSpectrumNumBasis;
        selectedModePairBatch.qSpectrumCondG = qSpectrumCondG;
        selectedModePairBatch.qSpectrumRankG = qSpectrumRankG;
        selectedModePairBatch.qSpectrumBasisInfo = qSpectrumBasisInfo;

        batchResults.selectedModePair = selectedModePairBatch;
        save_selected_mode_pair_batch_if_requested(selectedModePairBatch, cfg);
    end

    if saveStructureFit
        structureFitBatch = struct();
        structureFitBatch.A_trend_all = A_trend_all;
        structureFitBatch.A_oscillation_all = A_oscillation_all;
        structureFitBatch.A_trend_plus_oscillation_all = A_trend_plus_oscillation_all;
        structureFitBatch.q = q;
        structureFitBatch.t = t;
        structureFitBatch.trendModeIdx_all = trendModeIdx_all;
        structureFitBatch.oscillationModeIdx_all = oscillationModeIdx_all;
        structureFitBatch.trendRepresentativeMode = trendRepresentativeMode;
        structureFitBatch.oscillationRepresentativeMode = oscillationRepresentativeMode;
        structureFitBatch.trendPeriod_fs = trendPeriod_fs;
        structureFitBatch.trendDecayTime_fs = trendDecayTime_fs;
        structureFitBatch.oscillationPeriod_fs = oscillationPeriod_fs;
        structureFitBatch.oscillationDecayTime_fs = oscillationDecayTime_fs;
        structureFitBatch.oscillationWavenumber_cm1 = oscillationWavenumber_cm1;
        structureFitBatch.amplitudeMethod = 'global_delay_vector_least_squares';
        structureFitBatch.betaFitInfo = perRunBetaFitInfo;
        structureFitBatch.description = ...
            'De-Hankelized HODMD 2D components for structure fitting.';

        batchResults.structureFit = structureFitBatch;
        save_structure_fit_batch_if_requested(structureFitBatch, cfg);
    end

    fprintf('\nBatch target-mode summary:\n');
    disp(summaryTable);

    if isfield(batchResults.mode23Projection, 'summaryTable')
        fprintf('\nBatch mode(2,3) projection summary against zero-noise reference:\n');
        disp(batchResults.mode23Projection.summaryTable);
    end

    plot_batch_target_statistics(batchResults);
end


function singleResult = analyze_single_hodmd_target_metrics(A_slice, q, t, dt, cfg)
    A_work = preprocess_hodmd_matrix(A_slice, q, cfg);
    A_scale = ones(size(A_work, 1), 1);

    Uq = eye(size(A_work, 1));
    B = A_work;
    rState = size(B, 1);

    d = cfg.hodmd.delayOrder;

    if d == 0
        Z = B;
    else
        Z = make_delay_matrix(B, d);
    end

    if size(Z, 2) < 2
        error('DMD/HODMD requires at least two snapshots after preprocessing.');
    end

    Z1 = Z(:, 1:end-1);
    Z2 = Z(:, 2:end);

    [hodmd, dmdInfo] = exact_dmd(Z1, Z2, cfg.hodmd.dmdRank, dt); %#ok<ASGLU>
    [hodmd.beta, hodmd.b_firstSnapshot, hodmd.betaFitInfo] = ...
        fit_global_dmd_amplitudes(hodmd.Phi, hodmd.lambda, Z);
    hodmd.b = hodmd.beta;
    hodmd = sort_hodmd_modes_by_amplitude(hodmd);

    restrictedBetaComparison = struct();
    if isfield(cfg, 'restrictedBetaComparison') && ...
            cfg.restrictedBetaComparison.enabled
        restrictedBetaComparison = compare_restricted_global_beta(hodmd, Z);
    end

    projectionModePairA = [];
    projectionModePairIndices = [];
    if isfield(cfg, 'mode23Projection') && cfg.mode23Projection.enabled
        projectionModePairIndices = unique_valid_indices( ...
            cfg.mode23Projection.modeIdx, numel(hodmd.lambda));
        B_projectionPair = reconstruct_hodmd_modes( ...
            hodmd, projectionModePairIndices, rState, d, size(A_slice, 2));
        projectionModePairA = Uq * B_projectionPair;
        projectionModePairA = undo_component_scale( ...
            projectionModePairA, A_scale, cfg);
    end

    modeTable = build_hodmd_mode_table(hodmd, dt);
    dynamicModeTable = table();
    if isfield(cfg, 'batchTargetWavenumber_cm1')
        [targetGroups, dynamicModeTable] = select_nearest_dynamic_mode( ...
            hodmd, modeTable, cfg.batchTargetWavenumber_cm1);
    else
        targetGroups = find_target_hodmd_modes(hodmd, modeTable, cfg.target);
    end
    targetMetrics = extract_target_metrics(modeTable, targetGroups, cfg);
    if isfield(cfg, 'batchTargetWavenumber_cm1') && ~isempty(targetGroups{1})
        targetMetrics(1).amplitude = sum(modeTable.Amplitude(targetGroups{1}));
    end
    structureFitComponents = build_structure_fit_components( ...
        q, t, hodmd, Uq, A_scale, cfg, rState, d, size(A_slice, 2), modeTable);
    selectedModePair = build_selected_mode_pair_reconstruction( ...
        q, t, hodmd, Uq, A_scale, cfg, rState, d, size(A_slice, 2), modeTable);

    singleResult = struct();
    singleResult.modeTable = modeTable;
    singleResult.dynamicModeTable = dynamicModeTable;
    singleResult.targetGroups = targetGroups;
    singleResult.targetMetrics = targetMetrics;
    singleResult.structureFitComponents = structureFitComponents;
    singleResult.selectedModePair = selectedModePair;
    singleResult.betaFitInfo = hodmd.betaFitInfo;
    singleResult.restrictedBetaComparison = restrictedBetaComparison;
    singleResult.projectionModePairA = projectionModePairA;
    singleResult.projectionModePairIndices = projectionModePairIndices;
    singleResult.stateDimension = rState;
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


function [A_crop, q_crop, info] = crop_q_range(A_in, q_in, qRange)
    q_in = q_in(:);
    info = struct('applied', false, 'requestedRange', qRange, ...
        'selectedIndices', (1:numel(q_in)).');
    A_crop = A_in;
    q_crop = q_in;

    if isempty(qRange)
        return;
    end

    assert(isnumeric(qRange) && isreal(qRange) && numel(qRange) == 2 && ...
        all(isfinite(qRange(:))) && qRange(1) < qRange(2), ...
        'cfg.cropQRange must be [] or an increasing finite range [qMin qMax].');
    idxQ = find(q_in >= qRange(1) & q_in <= qRange(2));
    assert(numel(idxQ) >= 2, ...
        'cfg.cropQRange must retain at least two q points from the input.');

    A_crop = A_in(idxQ, :, :);
    q_crop = q_in(idxQ);
    info.applied = true;
    info.selectedIndices = idxQ;
end


function targetMetrics = extract_target_metrics(modeTable, targetGroups, cfg)
    nTargets = numel(targetGroups);

    targetMetrics = repmat(struct( ...
        'modeIndex', nan, ...
        'modeGroup', [], ...
        'modeGroupText', '[]', ...
        'frequency_fs', nan, ...
        'wavenumber_cm1', nan, ...
        'period_fs', nan, ...
        'growthRate_fs', nan, ...
        'decayTime_fs', nan, ...
        'amplitude', nan), nTargets, 1);

    for iTarget = 1:nTargets
        modeIdx = targetGroups{iTarget};
        modeIdx = unique_valid_indices(modeIdx, height(modeTable));

        if isempty(modeIdx)
            continue;
        end

        repMode = choose_representative_mode(modeIdx, modeTable);

        targetMetrics(iTarget).modeIndex = repMode;
        targetMetrics(iTarget).modeGroup = modeIdx;
        targetMetrics(iTarget).modeGroupText = mat2str(modeIdx);
        targetMetrics(iTarget).frequency_fs = modeTable.Frequency(repMode);
        targetMetrics(iTarget).wavenumber_cm1 = modeTable.Wavenumber_cm1(repMode);
        targetMetrics(iTarget).period_fs = modeTable.Period_fs(repMode);
        targetMetrics(iTarget).growthRate_fs = modeTable.GrowthRate(repMode);
        targetMetrics(iTarget).decayTime_fs = modeTable.DecayTime_fs(repMode);
        targetMetrics(iTarget).amplitude = modeTable.Amplitude(repMode);
    end
end


function summaryTable = build_batch_target_summary_table(targetPeriods, ...
    Frequency_fs, Wavenumber_cm1, Period_fs, GrowthRate_fs, DecayTime_fs, Amplitude)
    [FrequencyMean_fs, FrequencyStd_fs, NValid] = nan_column_stats(Frequency_fs);
    [WavenumberMean_cm1, WavenumberStd_cm1] = nan_column_stats(Wavenumber_cm1);
    [PeriodMean_fs, PeriodStd_fs] = nan_column_stats(Period_fs);
    [GrowthRateMean_fs, GrowthRateStd_fs] = nan_column_stats(GrowthRate_fs);
    [DecayTimeMean_fs, DecayTimeStd_fs] = nan_column_stats(DecayTime_fs);
    [AmplitudeMean, AmplitudeStd] = nan_column_stats(Amplitude);

    TargetPeriod_fs = targetPeriods(:);

    summaryTable = table(TargetPeriod_fs, NValid, ...
        FrequencyMean_fs, FrequencyStd_fs, ...
        WavenumberMean_cm1, WavenumberStd_cm1, ...
        PeriodMean_fs, PeriodStd_fs, ...
        GrowthRateMean_fs, GrowthRateStd_fs, ...
        DecayTimeMean_fs, DecayTimeStd_fs, ...
        AmplitudeMean, AmplitudeStd);
end


function [mu, sigma, nValid] = nan_column_stats(X)
    nCols = size(X, 2);
    mu = nan(nCols, 1);
    sigma = nan(nCols, 1);
    nValid = zeros(nCols, 1);

    for iCol = 1:nCols
        vals = X(:, iCol);
        vals = vals(isfinite(vals));
        nValid(iCol) = numel(vals);

        if isempty(vals)
            continue;
        end

        mu(iCol) = mean(vals);

        if isscalar(vals)
            sigma(iCol) = 0;
        else
            sigma(iCol) = std(vals, 0);
        end
    end
end


function truth = batch_truth_values(cfg)
    truth = struct();
    truth.enabled = false;
    truth.period_fs = nan;
    truth.wavenumber_cm1 = nan;
    truth.dampingTime_fs = nan;

    if ~isfield(cfg, 'mockI2Distance') || ~cfg.mockI2Distance.enabled
        return;
    end

    period_fs = cfg.mockI2Distance.period_fs;
    dampingTime_fs = cfg.mockI2Distance.tau_fs;

    if ~isfinite(period_fs) || period_fs <= 0
        return;
    end

    c_cm_per_fs = 2.99792458e-5;
    wavenumber_cm1 = (1 / period_fs) / c_cm_per_fs;

    truth.enabled = true;
    truth.period_fs = period_fs;
    truth.wavenumber_cm1 = wavenumber_cm1;
    truth.dampingTime_fs = dampingTime_fs;
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


function comparison = compare_restricted_global_beta(hodmd, Z)
    % Refit beta with selected mode columns only. Mode numbers are assumed
    % to follow the all-mode global-beta amplitude ordering.

    nModes = numel(hodmd.lambda);

    if nModes < 1
        error('Restricted beta comparison requires at least one DMD mode.');
    end

    mode1Idx = 1;
    mode123Idx = 1:min(3, nModes);
    allIdx = 1:nModes;

    [betaAll, fitAll] = fit_global_beta_for_mode_subset( ...
        hodmd.Phi, hodmd.lambda, Z, allIdx);
    [betaMode1, fitMode1] = fit_global_beta_for_mode_subset( ...
        hodmd.Phi, hodmd.lambda, Z, mode1Idx);
    [betaMode123, fitMode123] = fit_global_beta_for_mode_subset( ...
        hodmd.Phi, hodmd.lambda, Z, mode123Idx);

    betaMode1Expanded = complex(nan(nModes, 1), nan(nModes, 1));
    betaMode123Expanded = complex(nan(nModes, 1), nan(nModes, 1));
    betaMode1Expanded(mode1Idx) = betaMode1;
    betaMode123Expanded(mode123Idx) = betaMode123;

    coefficientTable = build_restricted_beta_coefficient_table( ...
        betaAll, betaMode1Expanded, betaMode123Expanded, mode123Idx);

    FitSet = ["all_modes"; "mode_1_only"; "modes_1_2_3"];
    NumModes = [numel(allIdx); numel(mode1Idx); numel(mode123Idx)];
    RelativeResidual = [fitAll.relErr; fitMode1.relErr; fitMode123.relErr];
    ResidualNorm = [fitAll.residualNorm; fitMode1.residualNorm; ...
        fitMode123.residualNorm];
    DesignRank = [fitAll.rankX; fitMode1.rankX; fitMode123.rankX];
    DesignCondition = [fitAll.condX; fitMode1.condX; fitMode123.condX];
    BetaNorm = [fitAll.betaNorm; fitMode1.betaNorm; fitMode123.betaNorm];

    fitQualityTable = table(FitSet, NumModes, RelativeResidual, ...
        ResidualNorm, DesignRank, DesignCondition, BetaNorm);

    comparison = struct();
    comparison.modeOrdering = ...
        'Modes sorted by modal amplitude from the all-mode global beta fit.';
    comparison.allModes.indices = allIdx;
    comparison.allModes.beta = betaAll;
    comparison.allModes.fitInfo = fitAll;
    comparison.mode1Only.indices = mode1Idx;
    comparison.mode1Only.beta = betaMode1;
    comparison.mode1Only.betaExpanded = betaMode1Expanded;
    comparison.mode1Only.fitInfo = fitMode1;
    comparison.modes123.indices = mode123Idx;
    comparison.modes123.beta = betaMode123;
    comparison.modes123.betaExpanded = betaMode123Expanded;
    comparison.modes123.fitInfo = fitMode123;
    comparison.coefficientTable = coefficientTable;
    comparison.fitQualityTable = fitQualityTable;
    comparison.allModeBetaConsistencyError = norm(betaAll - hodmd.beta(:));
end


function [beta, fitInfo] = fit_global_beta_for_mode_subset( ...
    Phi, lambda, Z, modeIdx)
    % Solve the global trajectory fit using only the requested DMD modes:
    %   vec(Z) ~= X_subset * beta_subset.

    nModes = numel(lambda);
    modeIdx = unique_valid_indices(modeIdx, nModes);

    if isempty(modeIdx)
        error('At least one valid mode is required for restricted beta fitting.');
    end

    lambdaSelected = lambda(modeIdx);
    PhiSelected = Phi(:, modeIdx);
    K = size(Z, 2);
    nSelected = numel(modeIdx);
    X = zeros(numel(Z), nSelected);

    for j = 1:nSelected
        temporalTrace = lambdaSelected(j).^(0:K-1);
        X(:, j) = reshape(PhiSelected(:, j) * temporalTrace, [], 1);
    end

    zVec = Z(:);
    beta = X \ zVec;
    residual = zVec - X * beta;
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

    fitInfo = struct();
    fitInfo.method = 'restricted_global_delay_vector_least_squares';
    fitInfo.modeIndices = modeIdx;
    fitInfo.numModes = nSelected;
    fitInfo.numSnapshots = K;
    fitInfo.numEquations = numel(Z);
    fitInfo.relErr = norm(residual) / max(norm(zVec), eps);
    fitInfo.residualNorm = norm(residual);
    fitInfo.rankX = rankX;
    fitInfo.condX = condX;
    fitInfo.betaNorm = norm(beta);
end


function T = build_restricted_beta_coefficient_table( ...
    betaAll, betaMode1Expanded, betaMode123Expanded, modeIdx)
    % Long-form table keeps complex beta values and their changes readable.

    betaSets = {betaAll, betaMode1Expanded, betaMode123Expanded};
    setNames = ["all_modes", "mode_1_only", "modes_1_2_3"];
    includeMask = {true(size(modeIdx)), modeIdx == 1, true(size(modeIdx))};
    nRows = sum(cellfun(@nnz, includeMask));

    FitSet = strings(nRows, 1);
    Mode = zeros(nRows, 1);
    BetaReal = zeros(nRows, 1);
    BetaImag = zeros(nRows, 1);
    BetaAbs = zeros(nRows, 1);
    BetaPhase_rad = zeros(nRows, 1);
    AbsChangeVsAll_pct = zeros(nRows, 1);
    PhaseChangeVsAll_rad = zeros(nRows, 1);

    iRow = 0;
    for iSet = 1:numel(betaSets)
        selectedRows = modeIdx(includeMask{iSet});

        for iMode = selectedRows
            iRow = iRow + 1;
            betaValue = betaSets{iSet}(iMode);
            betaReference = betaAll(iMode);

            FitSet(iRow) = setNames(iSet);
            Mode(iRow) = iMode;
            BetaReal(iRow) = real(betaValue);
            BetaImag(iRow) = imag(betaValue);
            BetaAbs(iRow) = abs(betaValue);
            BetaPhase_rad(iRow) = angle(betaValue);
            AbsChangeVsAll_pct(iRow) = 100 * ...
                (abs(betaValue) - abs(betaReference)) / ...
                max(abs(betaReference), eps);

            phaseDifference = angle(betaValue) - angle(betaReference);
            PhaseChangeVsAll_rad(iRow) = atan2( ...
                sin(phaseDifference), cos(phaseDifference));
        end
    end

    T = table(FitSet, Mode, BetaReal, BetaImag, BetaAbs, BetaPhase_rad, ...
        AbsChangeVsAll_pct, PhaseChangeVsAll_rad);
end


function print_restricted_global_beta_comparison(comparison)
    fprintf('\nRestricted global beta coefficient comparison:\n');
    fprintf('  Mode numbers follow the all-mode beta amplitude ordering.\n');
    fprintf('  Each restricted fit resolves its coefficients from all columns of Z.\n\n');
    disp(comparison.coefficientTable);

    fprintf('Restricted global beta fit quality over the full delay matrix Z:\n');
    disp(comparison.fitQualityTable);
    fprintf('  all-mode beta consistency error = %.6g\n\n', ...
        comparison.allModeBetaConsistencyError);
end


function reference = prepare_mode23_projection_reference(q, t, dt, cfg)
    reference = struct();
    reference.valid = false;
    reference.file = '';
    reference.A_mode23 = [];
    reference.modeIndices = [];
    reference.modeTable = table();

    referenceFile = resolve_mode23_projection_reference_file(cfg);
    reference.file = referenceFile;

    if isempty(referenceFile) || ~isfile(referenceFile)
        warning(['Mode(2,3) projection diagnostic skipped because the ', ...
            'zero-noise reference MAT file was not found: %s'], referenceFile);
        return;
    end

    [A_reference, q_reference, t_reference] = ...
        load_mode23_projection_reference_matrix(referenceFile, cfg);

    qTolerance = 1e-10 * max(1, max(abs(q)));
    tTolerance = 1e-10 * max(1, max(abs(t)));

    if numel(q_reference) ~= numel(q) || ...
            max(abs(q_reference(:) - q(:))) > qTolerance
        warning('Mode(2,3) projection skipped because reference q does not match input q.');
        return;
    end

    if numel(t_reference) ~= numel(t) || ...
            max(abs(t_reference(:) - t(:))) > tTolerance
        warning('Mode(2,3) projection skipped because reference t does not match input t.');
        return;
    end

    referenceCfg = cfg;
    referenceCfg.restrictedBetaComparison.enabled = false;
    singleReference = analyze_single_hodmd_target_metrics( ...
        A_reference, q_reference, t_reference, dt, referenceCfg);

    if isempty(singleReference.projectionModePairA)
        warning('Mode(2,3) projection skipped because the reference pair was not reconstructed.');
        return;
    end

    reference.valid = true;
    reference.A_mode23 = real(singleReference.projectionModePairA);
    reference.modeIndices = singleReference.projectionModePairIndices;
    reference.modeTable = singleReference.modeTable;

    representativeMode = choose_representative_mode( ...
        reference.modeIndices, reference.modeTable);
    reference.representativeMode = representativeMode;
    reference.period_fs = reference.modeTable.Period_fs(representativeMode);
    reference.wavenumber_cm1 = ...
        reference.modeTable.Wavenumber_cm1(representativeMode);
    reference.decayTime_fs = ...
        reference.modeTable.DecayTime_fs(representativeMode);
end


function referenceFile = resolve_mode23_projection_reference_file(cfg)
    configuredFile = '';
    if isfield(cfg.mode23Projection, 'referenceMatFile')
        configuredFile = char(string(cfg.mode23Projection.referenceMatFile));
    end

    if ~isempty(strtrim(configuredFile))
        referenceFile = configuredFile;
        return;
    end

    [inputDir, inputName, inputExt] = fileparts(char(cfg.matFile));
    referenceName = regexprep(inputName, 'noise\d+', 'noise0', 'once');

    if strcmp(referenceName, inputName)
        referenceFile = '';
        warning(['Could not infer a zero-noise reference filename from ', ...
            'cfg.matFile. Set cfg.mode23Projection.referenceMatFile explicitly.']);
        return;
    end

    referenceFile = fullfile(inputDir, [referenceName inputExt]);
end


function [A_reference, q_reference, t_reference_fs] = ...
        load_mode23_projection_reference_matrix(referenceFile, cfg)
    S = load(referenceFile);

    assert(isfield(S, cfg.varA) && isfield(S, cfg.varq) && ...
        isfield(S, cfg.vart), ...
        'Reference MAT file must contain %s, %s, and %s.', ...
        cfg.varA, cfg.varq, cfg.vart);

    A_reference = double(S.(cfg.varA));
    q_reference = double(S.(cfg.varq));
    t_reference = double(S.(cfg.vart));
    q_reference = q_reference(:);
    t_reference = t_reference(:).';

    if size(A_reference, 1) == numel(q_reference) && ...
            size(A_reference, 2) == numel(t_reference)
        % Already q x time, with an optional series dimension.
    elseif size(A_reference, 1) == numel(t_reference) && ...
            size(A_reference, 2) == numel(q_reference)
        if ismatrix(A_reference)
            A_reference = A_reference.';
        else
            A_reference = permute(A_reference, [2 1 3]);
        end
    else
        error('Cannot orient zero-noise reference A to q x time.');
    end

    if ndims(A_reference) == 3
        A_reference = A_reference(:, :, 1);
    end

    [A_reference, q_reference] = crop_q_range( ...
        A_reference, q_reference, cfg.cropQRange);

    [t_reference, sortIdx] = sort(t_reference, 'ascend');
    A_reference = A_reference(:, sortIdx);

    idx0 = cfg.cropStartIndex;
    assert(idx0 >= 1 && idx0 <= numel(t_reference), ...
        'cfg.cropStartIndex is invalid for the zero-noise reference.');
    t_reference = t_reference(idx0:end);
    A_reference = A_reference(:, idx0:end);
    t_reference_fs = t_reference .* cfg.fsPerInputTimeUnit;


    assert(all(isfinite(A_reference(:))), ...
        'Zero-noise reference A contains Inf or unreplaced NaN values.');
end


function diagnostics = compute_mode23_projection_diagnostics( ...
        A_mode23, A_referenceMode23, modeIndices, reference)
    target = real(A_mode23);
    template = real(A_referenceMode23);
    valid = isfinite(target) & isfinite(template);

    assert(any(valid(:)), ...
        'No finite points are available for mode(2,3) projection.');

    targetValues = target(valid);
    templateValues = template(valid);
    templateEnergy = sum(templateValues .^ 2);
    assert(templateEnergy > eps, ...
        'The zero-noise mode(2,3) reference has zero numerical energy.');

    projectionCoefficient = ...
        sum(templateValues .* targetValues) ./ templateEnergy;
    parallelValues = projectionCoefficient .* templateValues;
    orthogonalValues = targetValues - parallelValues;

    targetNorm = norm(targetValues);
    templateNorm = sqrt(templateEnergy);
    parallelNorm = norm(parallelValues);
    orthogonalNorm = norm(orthogonalValues);
    targetNormSafe = max(targetNorm, eps);

    A_parallel = nan(size(target));
    A_orthogonal = nan(size(target));
    A_parallel(valid) = parallelValues;
    A_orthogonal(valid) = orthogonalValues;

    diagnostics = struct();
    diagnostics.modeIndices = modeIndices;
    diagnostics.referenceFile = reference.file;
    diagnostics.referenceModeIndices = reference.modeIndices;
    diagnostics.referencePeriod_fs = reference.period_fs;
    diagnostics.referenceWavenumber_cm1 = reference.wavenumber_cm1;
    diagnostics.referenceDecayTime_fs = reference.decayTime_fs;
    diagnostics.projectionCoefficient = projectionCoefficient;
    diagnostics.targetNorm = targetNorm;
    diagnostics.referenceNorm = templateNorm;
    diagnostics.parallelNorm = parallelNorm;
    diagnostics.orthogonalNorm = orthogonalNorm;
    diagnostics.targetToReferenceNormRatio = targetNorm ./ templateNorm;
    diagnostics.parallelFraction = parallelNorm ./ targetNormSafe;
    diagnostics.orthogonalFraction = orthogonalNorm ./ targetNormSafe;
    diagnostics.parallelEnergyFraction = diagnostics.parallelFraction .^ 2;
    diagnostics.orthogonalEnergyFraction = diagnostics.orthogonalFraction .^ 2;
    diagnostics.cosineSimilarity = sum(templateValues .* targetValues) ./ ...
        max(templateNorm .* targetNorm, eps);
    diagnostics.orthogonalityCheck = sum(parallelValues .* orthogonalValues) ./ ...
        max(parallelNorm .* orthogonalNorm, eps);
    diagnostics.A_mode23 = target;
    diagnostics.A_referenceMode23 = template;
    diagnostics.A_parallel = A_parallel;
    diagnostics.A_orthogonal = A_orthogonal;
    diagnostics.validMask = valid;
end


function print_mode23_projection_diagnostics(diagnostics)
    fprintf('\nMode(2,3) projection against zero-noise HODMD reference:\n');
    fprintf('  current modes                     = %s\n', ...
        mat2str(diagnostics.modeIndices));
    fprintf('  reference modes                   = %s\n', ...
        mat2str(diagnostics.referenceModeIndices));
    fprintf('  reference period / tau            = %.6g / %.6g fs\n', ...
        diagnostics.referencePeriod_fs, diagnostics.referenceDecayTime_fs);
    fprintf('  reference wavenumber              = %.6g cm^-1\n', ...
        diagnostics.referenceWavenumber_cm1);
    fprintf('  ||A_mode23||_F                    = %.12g\n', ...
        diagnostics.targetNorm);
    fprintf('  ||A_mode23_zeroNoise||_F          = %.12g\n', ...
        diagnostics.referenceNorm);
    fprintf('  projection coefficient            = %.12g\n', ...
        diagnostics.projectionCoefficient);
    fprintf('  ||A_parallel||_F                  = %.12g\n', ...
        diagnostics.parallelNorm);
    fprintf('  ||A_orthogonal||_F                = %.12g\n', ...
        diagnostics.orthogonalNorm);
    fprintf('  parallel energy fraction          = %.6f\n', ...
        diagnostics.parallelEnergyFraction);
    fprintf('  orthogonal energy fraction        = %.6f\n', ...
        diagnostics.orthogonalEnergyFraction);
    fprintf('  target / zero-noise norm          = %.6f\n', ...
        diagnostics.targetToReferenceNormRatio);
    fprintf('  cosine similarity                 = %.6f\n', ...
        diagnostics.cosineSimilarity);
    fprintf('  reference file                    = %s\n\n', ...
        diagnostics.referenceFile);
end


function summaryTable = build_mode23_projection_batch_summary(batch)
    Metric = ["Mode23Norm"; "ParallelNorm"; "OrthogonalNorm"; ...
        "ProjectionCoefficient"; "ParallelEnergyFraction"; ...
        "OrthogonalEnergyFraction"; "CosineSimilarity"];
    values = {batch.targetNorm, batch.parallelNorm, batch.orthogonalNorm, ...
        batch.projectionCoefficient, batch.parallelEnergyFraction, ...
        batch.orthogonalEnergyFraction, batch.cosineSimilarity};

    Mean = nan(numel(values), 1);
    Std = nan(numel(values), 1);
    N = zeros(numel(values), 1);

    for i = 1:numel(values)
        finiteValues = values{i}(isfinite(values{i}));
        N(i) = numel(finiteValues);
        if ~isempty(finiteValues)
            Mean(i) = mean(finiteValues);
            Std(i) = std(finiteValues, 0);
        end
    end

    summaryTable = table(Metric, Mean, Std, N);
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


function modeComponentRecon = reconstruct_all_mode_components( ...
    hodmd, Uq, A_scale, cfg, r, d, Nt, modeTable)
    groups = group_modes_by_conjugacy(hodmd.lambda);
    % Rank AFTER conjugate pairing. Keep eigencomponent indices unchanged.
    groupAmplitude = arrayfun(@(g) sum(modeTable.Amplitude(g.modeIdx)), groups);
    [~, order] = sort(groupAmplitude, 'descend');
    groups = groups(order);
    nComp = numel(groups);
    nQ = size(Uq, 1);

    components = repmat(struct( ...
        'componentIndex', [], ...
        'dynamicMode', '', ...
        'type', '', ...
        'isComplex', false, ...
        'modeIdx', [], ...
        'representativeMode', [], ...
        'frequency_fs', [], ...
        'wavenumber_cm1', [], ...
        'period_fs', [], ...
        'decayTime_fs', [], ...
        'amplitude', [], ...
        'A', []), nComp, 1);

    A_modeComponents = zeros(nQ, Nt, nComp);

    Component = (1:nComp).';
    DynamicMode = strings(nComp, 1);
    dynamicRank = 0;
    Type = cell(nComp, 1);
    ModeIndices = cell(nComp, 1);
    RepresentativeMode = zeros(nComp, 1);
    Frequency = zeros(nComp, 1);
    Wavenumber_cm1 = zeros(nComp, 1);
    Period_fs = zeros(nComp, 1);
    DecayTime_fs = zeros(nComp, 1);
    Amplitude = zeros(nComp, 1);

    for k = 1:nComp
        modeIdx = groups(k).modeIdx;

        B_component = reconstruct_hodmd_modes(hodmd, modeIdx, r, d, Nt);
        A_component = Uq * B_component;
        A_component = undo_component_scale(A_component, A_scale, cfg);

        repMode = choose_representative_mode(modeIdx, modeTable);

        components(k).componentIndex = k;
        if strcmp(groups(k).type, 'complex_unpaired')
            DynamicMode(k) = "Unpaired";
        else
            dynamicRank = dynamicRank + 1;
            DynamicMode(k) = sprintf('M%d', dynamicRank);
        end
        components(k).dynamicMode = char(DynamicMode(k));
        components(k).type = groups(k).type;
        components(k).isComplex = groups(k).isComplex;
        components(k).modeIdx = modeIdx;
        components(k).representativeMode = repMode;
        components(k).frequency_fs = modeTable.Frequency(repMode);
        components(k).wavenumber_cm1 = modeTable.Wavenumber_cm1(repMode);
        components(k).period_fs = modeTable.Period_fs(repMode);
        components(k).decayTime_fs = modeTable.DecayTime_fs(repMode);
        components(k).amplitude = sum(modeTable.Amplitude(modeIdx));
        components(k).A = A_component;

        A_modeComponents(:, :, k) = A_component;

        Type{k} = groups(k).type;
        ModeIndices{k} = mat2str(modeIdx);
        RepresentativeMode(k) = repMode;
        Frequency(k) = components(k).frequency_fs;
        Wavenumber_cm1(k) = components(k).wavenumber_cm1;
        Period_fs(k) = components(k).period_fs;
        DecayTime_fs(k) = components(k).decayTime_fs;
        Amplitude(k) = components(k).amplitude;
    end

    componentTable = table(Component, DynamicMode, Type, ModeIndices, ...
        RepresentativeMode, Frequency, Wavenumber_cm1, Period_fs, ...
        DecayTime_fs, Amplitude);

    modeComponentRecon = struct();
    modeComponentRecon.components = components;
    modeComponentRecon.A = A_modeComponents;
    modeComponentRecon.table = componentTable;
end


function [targetGroups, dynamicModeTable] = select_nearest_dynamic_mode(hodmd, modeTable, target_cm1)
    % An eigencomponent belongs to one eigenvalue. A dynamic mode is either
    % one real eigencomponent or one complete complex-conjugate pair.
    groups = group_modes_by_conjugacy(hodmd.lambda);
    complete = ~strcmp({groups.type}, 'complex_unpaired');
    if any(~complete)
        warning('Unpaired complex eigencomponents excluded from dynamic-mode ranking.');
    end
    groups = groups(complete);
    n = numel(groups);
    EigencomponentIndices = cell(n, 1);
    Type = cell(n, 1);
    Amplitude = zeros(n, 1);
    Wavenumber_cm1 = nan(n, 1);
    for k = 1:n
        EigencomponentIndices{k} = groups(k).modeIdx;
        Type{k} = groups(k).type;
        Amplitude(k) = sum(modeTable.Amplitude(groups(k).modeIdx));
        rep = choose_representative_mode(groups(k).modeIdx, modeTable);
        Wavenumber_cm1(k) = abs(modeTable.Wavenumber_cm1(rep));
    end
    [Amplitude, order] = sort(Amplitude, 'descend');
    EigencomponentIndices = EigencomponentIndices(order);
    Type = Type(order);
    Wavenumber_cm1 = Wavenumber_cm1(order);
    DynamicMode = compose('M%d', (1:n).');
    DistanceToTarget_cm1 = abs(Wavenumber_cm1 - target_cm1);
    Selected = false(n, 1);
    targetGroups = {[]};
    candidates = find(strcmp(Type, 'complex_pair') & isfinite(Wavenumber_cm1));
    if ~isempty(candidates)
        [~, nearest] = min(DistanceToTarget_cm1(candidates));
        selected = candidates(nearest);
        Selected(selected) = true;
        targetGroups{1} = EigencomponentIndices{selected};
    end
    dynamicModeTable = table(DynamicMode, EigencomponentIndices, Type, ...
        Amplitude, Wavenumber_cm1, DistanceToTarget_cm1, Selected);
end


function groups = group_modes_by_conjugacy(lambda)
    lambda = lambda(:);
    nModes = numel(lambda);
    used = false(nModes, 1);
    tol = 1e-6;

    groups = repmat(struct( ...
        'modeIdx', [], ...
        'type', '', ...
        'isComplex', false), 0, 1);

    for k = 1:nModes
        if used(k)
            continue;
        end

        pairIdx = find_conjugate_partner(lambda, k, tol);

        if isempty(pairIdx)
            modeIdx = k;
            if abs(imag(lambda(k))) <= tol * max(1, abs(lambda(k)))
                type = 'real';
                isComplex = false;
            else
                type = 'complex_unpaired';
                isComplex = true;
            end
        else
            modeIdx = unique([k pairIdx], 'stable');
            type = 'complex_pair';
            isComplex = true;
        end

        used(modeIdx) = true;

        groups(end+1, 1).modeIdx = modeIdx; %#ok<AGROW>
        groups(end).type = type;
        groups(end).isComplex = isComplex;
    end
end


function repMode = choose_representative_mode(modeIdx, modeTable)
    idx = modeIdx(:).';

    if isscalar(idx)
        repMode = idx;
        return;
    end

    posFreq = idx(modeTable.Frequency(idx) > 0);

    if ~isempty(posFreq)
        [~, best] = max(modeTable.Amplitude(posFreq));
        repMode = posFreq(best);
    else
        [~, best] = max(modeTable.Amplitude(idx));
        repMode = idx(best);
    end
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


function betaComparison = compare_b_beta_mode_pair( ...
        hodmd_b, hodmd_beta, q, t, Uq, A_scale, cfg, r, d, Nt)
    modeIdx = cfg.compareBeta.modeIdx(:).';
    modeIdx = unique_valid_indices(modeIdx, numel(hodmd_beta.lambda));

    if isempty(modeIdx)
        warning('b_j vs beta_j comparison skipped because cfg.compareBeta.modeIdx is empty or invalid.');
        betaComparison = struct();
        return;
    end

    Y_b = reconstruct_hodmd_modes(hodmd_b, modeIdx, r, d, Nt);
    A_b = Uq * Y_b;
    A_b = undo_component_scale(A_b, A_scale, cfg);

    Y_beta = reconstruct_hodmd_modes(hodmd_beta, modeIdx, r, d, Nt);
    A_beta = Uq * Y_beta;
    A_beta = undo_component_scale(A_beta, A_scale, cfg);

    A_diff = A_beta - A_b;

    normB = norm(A_b, 'fro');
    normBeta = norm(A_beta, 'fro');
    normDiff = norm(A_diff, 'fro');
    denomCorr = max(normB * normBeta, eps);
    corrComplex = sum(conj(A_b(:)) .* A_beta(:)) / denomCorr;

    bSel = hodmd_b.b(modeIdx);
    betaSel = hodmd_beta.b(modeIdx);
    absRatio = abs(betaSel) ./ max(abs(bSel), eps);

    phaseDiff = nan(size(modeIdx(:)));
    nonzeroB = abs(bSel(:)) > eps;
    phaseDiff(nonzeroB) = angle(betaSel(nonzeroB) ./ bSel(nonzeroB));

    Mode = modeIdx(:);
    Real_bj = real(bSel(:));
    Imag_bj = imag(bSel(:));
    Abs_bj = abs(bSel(:));
    Real_betaj = real(betaSel(:));
    Imag_betaj = imag(betaSel(:));
    Abs_betaj = abs(betaSel(:));
    AbsRatio_beta_over_b = absRatio(:);
    PhaseDiff_rad = phaseDiff(:);

    coefficientTable = table(Mode, ...
        Real_bj, Imag_bj, Abs_bj, ...
        Real_betaj, Imag_betaj, Abs_betaj, ...
        AbsRatio_beta_over_b, PhaseDiff_rad);

    fitInfo = hodmd_beta.betaFitInfo;

    betaComparison = struct();
    betaComparison.modeIdx = modeIdx;
    betaComparison.A_b = A_b;
    betaComparison.A_beta = A_beta;
    betaComparison.A_diff = A_diff;
    betaComparison.Y_b = Y_b;
    betaComparison.Y_beta = Y_beta;
    betaComparison.coefficientTable = coefficientTable;
    betaComparison.relDiff_beta_minus_b_over_b = normDiff / max(normB, eps);
    betaComparison.relDiff_beta_minus_b_over_beta = normDiff / max(normBeta, eps);
    betaComparison.correlation = corrComplex;
    betaComparison.fitInfo = fitInfo;

    fprintf('\nb_j vs beta_j comparison for mode pair %s\n', mat2str(modeIdx));
    fprintf('  first-delay-vector b_j rel. error over Z = %.6g\n', ...
        fitInfo.relErr_firstSnapshotB);
    fprintf('  global beta_j rel. error over Z          = %.6g\n', ...
        fitInfo.relErr_globalBeta);
    fprintf('  ||A_beta - A_b||_F / ||A_b||_F          = %.6g\n', ...
        betaComparison.relDiff_beta_minus_b_over_b);
    fprintf('  ||A_beta - A_b||_F / ||A_beta||_F       = %.6g\n', ...
        betaComparison.relDiff_beta_minus_b_over_beta);
    fprintf('  2D reconstruction correlation           = %.6g%+.6gi\n\n', ...
        real(corrComplex), imag(corrComplex));
    disp(coefficientTable);

    if isfield(cfg.compareBeta, 'plotResults') && cfg.compareBeta.plotResults
        plot_b_beta_mode_pair_comparison(t, q, A_b, A_beta, A_diff, ...
            modeIdx, betaComparison, cfg);
        plot_b_beta_coefficient_comparison(modeIdx, bSel, betaSel, cfg);
    end
end


function plot_b_beta_mode_pair_comparison(t, q, A_b, A_beta, A_diff, ...
        modeIdx, betaComparison, cfg)
    figName = sprintf('b_j vs beta_j mode %s 2D reconstruction', ...
        mode_index_suffix(modeIdx));
    figure('Color','w','Name',figName, 'WindowStyle','docked');

    tiledlayout(1, 3, 'TileSpacing','compact', 'Padding','compact');

    mainClim = compute_qt_clim([A_b(:); A_beta(:)], cfg);
    diffClim = compute_qt_clim(A_diff(:), cfg);

    nexttile;
    imagesc(q, t, A_b.');
    set(gca, 'YDir', 'reverse');
    axis square;
    xlabel('q');
    ylabel('time [fs]');
    title('b_j: first delay vector');
    colorbar;
    if ~isempty(mainClim)
        clim(mainClim);
    end

    nexttile;
    imagesc(q, t, A_beta.');
    set(gca, 'YDir', 'reverse');
    axis square;
    xlabel('q');
    ylabel('time [fs]');
    title('beta_j: global fit');
    colorbar;
    if ~isempty(mainClim)
        clim(mainClim);
    end

    nexttile;
    imagesc(q, t, A_diff.');
    set(gca, 'YDir', 'reverse');
    axis square;
    xlabel('q');
    ylabel('time [fs]');
    title('beta_j - b_j');
    colorbar;
    if ~isempty(diffClim)
        clim(diffClim);
    end

    colormap(turbo);

    sgtitle(sprintf(['mode %s pair: rel diff over beta = %.3g, ' ...
        'corr = %.4g%+.4gi'], ...
        mat2str(modeIdx), ...
        betaComparison.relDiff_beta_minus_b_over_beta, ...
        real(betaComparison.correlation), ...
        imag(betaComparison.correlation)), ...
        'Interpreter','none');
end


function plot_b_beta_coefficient_comparison(modeIdx, bSel, betaSel, ~)
    figure('Color','w','Name','b_j vs beta_j coefficients', ...
        'WindowStyle','docked');

    tiledlayout(1, 3, 'TileSpacing','compact', 'Padding','compact');

    labels = categorical(compose('mode %d', modeIdx));
    labels = reordercats(labels, cellstr(compose('mode %d', modeIdx)));

    nexttile;
    bar(labels, [abs(bSel(:)), abs(betaSel(:))]);
    ylabel('coefficient magnitude');
    title('|b_j| vs |beta_j|');
    legend({'b_j','beta_j'}, 'Location','best');
    grid on;

    phaseDiff = nan(size(modeIdx(:)));
    nonzeroB = abs(bSel(:)) > eps;
    phaseDiff(nonzeroB) = angle(betaSel(nonzeroB) ./ bSel(nonzeroB));

    nexttile;
    stem(modeIdx, phaseDiff, 'filled', 'LineWidth', 1.2);
    xlabel('mode index');
    ylabel('phase(beta_j / b_j) [rad]');
    title('coefficient phase difference');
    grid on;

    nexttile;
    plot(real(bSel), imag(bSel), 'bo', 'MarkerFaceColor','b', ...
        'MarkerSize', 7);
    hold on;
    plot(real(betaSel), imag(betaSel), 'rs', 'MarkerFaceColor','r', ...
        'MarkerSize', 7);
    for k = 1:numel(modeIdx)
        plot([real(bSel(k)), real(betaSel(k))], ...
             [imag(bSel(k)), imag(betaSel(k))], 'k-', ...
             'LineWidth', 0.8);
        text(real(betaSel(k)), imag(betaSel(k)), ...
            sprintf(' %d', modeIdx(k)), 'Interpreter','none');
    end
    hold off;
    axis equal;
    grid on;
    xlabel('real coefficient');
    ylabel('imag coefficient');
    title('complex coefficient plane');
    legend({'b_j','beta_j'}, 'Location','best');
end


function Aout = undo_component_scale(Ain, ~, ~)
    % No q-wise normalization is applied in this publication workflow.
    Aout = Ain;
end


function tf = is_structure_fit_save_enabled(cfg)
    tf = isfield(cfg, 'structureFitSave') && ...
        isfield(cfg.structureFitSave, 'enabled') && ...
        cfg.structureFitSave.enabled;
end


function components = build_structure_fit_components( ...
    q, t, hodmd, Uq, A_scale, cfg, r, d, Nt, modeTable)
    components = struct();
    components.enabled = false;
    components.valid = false;
    components.q = q(:);
    components.t = t(:).';
    components.A_trend = [];
    components.A_oscillation = [];
    components.A_trend_plus_oscillation = [];
    components.trendModeIdx = [];
    components.oscillationModeIdx = [];
    components.trendRepresentativeMode = [];
    components.oscillationRepresentativeMode = [];
    components.trendPeriod_fs = nan;
    components.trendDecayTime_fs = nan;
    components.oscillationPeriod_fs = nan;
    components.oscillationDecayTime_fs = nan;
    components.oscillationWavenumber_cm1 = nan;
    components.trendModeTable = table();
    components.oscillationModeTable = table();
    components.amplitudeMethod = 'global_delay_vector_least_squares';
    components.betaFitInfo = struct();
    components.description = ...
        'De-Hankelized HODMD 2D components for structure fitting.';

    if ~is_structure_fit_save_enabled(cfg)
        return;
    end

    components.enabled = true;

    trendModeIdx = cfg.structureFitSave.trendModeIdx;
    oscillationModeIdx = cfg.structureFitSave.oscillationModeIdx;

    trendModeIdx = unique_valid_indices(trendModeIdx, numel(hodmd.lambda));
    oscillationModeIdx = unique_valid_indices(oscillationModeIdx, numel(hodmd.lambda));

    if isempty(trendModeIdx) && isempty(oscillationModeIdx)
        warning('Structure-fit component save skipped because all requested mode indices are invalid.');
        return;
    end

    B_trend = reconstruct_hodmd_modes(hodmd, trendModeIdx, r, d, Nt);
    A_trend = Uq * B_trend;
    A_trend = undo_component_scale(A_trend, A_scale, cfg);

    B_oscillation = reconstruct_hodmd_modes(hodmd, oscillationModeIdx, r, d, Nt);
    A_oscillation = Uq * B_oscillation;
    A_oscillation = undo_component_scale(A_oscillation, A_scale, cfg);

    components.valid = true;
    components.A_trend = A_trend;
    components.A_oscillation = A_oscillation;
    components.A_trend_plus_oscillation = A_trend + A_oscillation;
    components.betaFitInfo = hodmd.betaFitInfo;
    components.trendModeIdx = trendModeIdx;
    components.oscillationModeIdx = oscillationModeIdx;

    if ~isempty(trendModeIdx)
        repTrend = choose_representative_mode(trendModeIdx, modeTable);
        components.trendRepresentativeMode = repTrend;
        components.trendPeriod_fs = modeTable.Period_fs(repTrend);
        components.trendDecayTime_fs = modeTable.DecayTime_fs(repTrend);
        components.trendModeTable = modeTable(trendModeIdx, :);
    end

    if ~isempty(oscillationModeIdx)
        repOsc = choose_representative_mode(oscillationModeIdx, modeTable);
        components.oscillationRepresentativeMode = repOsc;
        components.oscillationPeriod_fs = modeTable.Period_fs(repOsc);
        components.oscillationDecayTime_fs = modeTable.DecayTime_fs(repOsc);
        components.oscillationWavenumber_cm1 = modeTable.Wavenumber_cm1(repOsc);
        components.oscillationModeTable = modeTable(oscillationModeIdx, :);
    end
end


function saveInfo = save_structure_fit_components_if_requested( ...
    q, t, hodmd, Uq, A_scale, cfg, r, d, Nt, modeTable)
    saveInfo = struct();
    saveInfo.enabled = false;
    saveInfo.saved = false;
    saveInfo.file = '';

    if ~is_structure_fit_save_enabled(cfg)
        return;
    end

    saveInfo.enabled = true;

    components = build_structure_fit_components( ...
        q, t, hodmd, Uq, A_scale, cfg, r, d, Nt, modeTable);

    if ~components.valid
        warning('Structure-fit component save skipped because requested modes were invalid.');
        return;
    end

    outFile = structure_fit_output_file(cfg);

    A_trend = components.A_trend;
    A_oscillation = components.A_oscillation;
    A_trend_plus_oscillation = components.A_trend_plus_oscillation;
    q = components.q;
    t = components.t;
    trendModeIdx = components.trendModeIdx;
    oscillationModeIdx = components.oscillationModeIdx;
    trendRepresentativeMode = components.trendRepresentativeMode;
    oscillationRepresentativeMode = components.oscillationRepresentativeMode;
    trendPeriod_fs = components.trendPeriod_fs;
    trendDecayTime_fs = components.trendDecayTime_fs;
    oscillationPeriod_fs = components.oscillationPeriod_fs;
    oscillationDecayTime_fs = components.oscillationDecayTime_fs;
    oscillationWavenumber_cm1 = components.oscillationWavenumber_cm1;
    trendModeTable = components.trendModeTable;
    oscillationModeTable = components.oscillationModeTable;
    amplitudeMethod = components.amplitudeMethod;
    betaFitInfo = components.betaFitInfo;
    description = components.description;

    save(outFile, ...
        'A_trend', 'A_oscillation', 'A_trend_plus_oscillation', ...
        'q', 't', ...
        'trendModeIdx', 'oscillationModeIdx', ...
        'trendRepresentativeMode', 'oscillationRepresentativeMode', ...
        'trendPeriod_fs', 'trendDecayTime_fs', ...
        'oscillationPeriod_fs', 'oscillationDecayTime_fs', ...
        'oscillationWavenumber_cm1', ...
        'trendModeTable', 'oscillationModeTable', ...
        'amplitudeMethod', 'betaFitInfo', ...
        'description');
    append_input_fwhm_if_available(outFile, cfg);

    saveInfo.saved = true;
    saveInfo.file = outFile;
    saveInfo.trendModeIdx = trendModeIdx;
    saveInfo.oscillationModeIdx = oscillationModeIdx;
    saveInfo.trendRepresentativeMode = trendRepresentativeMode;
    saveInfo.oscillationRepresentativeMode = oscillationRepresentativeMode;
    saveInfo.oscillationPeriod_fs = oscillationPeriod_fs;
    saveInfo.oscillationDecayTime_fs = oscillationDecayTime_fs;
    saveInfo.oscillationWavenumber_cm1 = oscillationWavenumber_cm1;

    fprintf('\nSaved structure-fit 2D HODMD components:\n');
    fprintf('  trend modes = %s\n', mat2str(trendModeIdx));
    fprintf('  oscillation modes = %s\n', mat2str(oscillationModeIdx));
    fprintf('  oscillation period = %.6g fs\n', oscillationPeriod_fs);
    fprintf('  oscillation damping time = %.6g fs\n', oscillationDecayTime_fs);
    fprintf('  file = %s\n\n', outFile);
end


function save_structure_fit_batch_if_requested(structureFitBatch, cfg)
    if ~is_structure_fit_save_enabled(cfg)
        return;
    end

    outFile = structure_fit_output_file(cfg);

    A_trend_all = structureFitBatch.A_trend_all;
    A_oscillation_all = structureFitBatch.A_oscillation_all;
    A_trend_plus_oscillation_all = structureFitBatch.A_trend_plus_oscillation_all;
    q = structureFitBatch.q;
    t = structureFitBatch.t;
    trendModeIdx_all = structureFitBatch.trendModeIdx_all;
    oscillationModeIdx_all = structureFitBatch.oscillationModeIdx_all;
    trendRepresentativeMode = structureFitBatch.trendRepresentativeMode;
    oscillationRepresentativeMode = structureFitBatch.oscillationRepresentativeMode;
    trendPeriod_fs = structureFitBatch.trendPeriod_fs;
    trendDecayTime_fs = structureFitBatch.trendDecayTime_fs;
    oscillationPeriod_fs = structureFitBatch.oscillationPeriod_fs;
    oscillationDecayTime_fs = structureFitBatch.oscillationDecayTime_fs;
    oscillationWavenumber_cm1 = structureFitBatch.oscillationWavenumber_cm1;
    amplitudeMethod = structureFitBatch.amplitudeMethod;
    betaFitInfo = structureFitBatch.betaFitInfo;
    description = structureFitBatch.description;

    save(outFile, ...
        'A_trend_all', 'A_oscillation_all', ...
        'A_trend_plus_oscillation_all', ...
        'q', 't', ...
        'trendModeIdx_all', 'oscillationModeIdx_all', ...
        'trendRepresentativeMode', 'oscillationRepresentativeMode', ...
        'trendPeriod_fs', 'trendDecayTime_fs', ...
        'oscillationPeriod_fs', 'oscillationDecayTime_fs', ...
        'oscillationWavenumber_cm1', ...
        'amplitudeMethod', 'betaFitInfo', ...
        'description', ...
        '-v7.3');
    append_input_fwhm_if_available(outFile, cfg);

    fprintf('\nSaved structure-fit 3D HODMD components:\n');
    fprintf('  trend modes = %s\n', mat2str(cfg.structureFitSave.trendModeIdx));
    fprintf('  oscillation modes = %s\n', mat2str(cfg.structureFitSave.oscillationModeIdx));
    fprintf('  file = %s\n\n', outFile);
end


function outFile = structure_fit_output_file(cfg)
    outputDir = global_beta_output_dir(cfg);

    if isfield(cfg.structureFitSave, 'outputDir') && ...
            ~isempty(cfg.structureFitSave.outputDir)
        outputDir = cfg.structureFitSave.outputDir;
    end

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    outputFile = input_based_output_file(cfg, 'globalBeta_structureFit');

    if isfield(cfg.structureFitSave, 'outputFile') && ...
            ~isempty(cfg.structureFitSave.outputFile)
        outputFile = configured_or_input_output_file( ...
            cfg, cfg.structureFitSave.outputFile, 'globalBeta_structureFit');
    end

    outFile = fullfile_if_no_folder(outputDir, outputFile);
end


function outputDir = global_beta_output_dir(cfg)
    outputDir = '';

    if isfield(cfg, 'outputDir') && ~isempty(cfg.outputDir)
        outputDir = cfg.outputDir;
    end

    if isempty(outputDir)
        scriptPath = mfilename('fullpath');
        outputDir = fileparts(scriptPath);

        if isempty(outputDir)
            outputDir = pwd;
        end

        outputDir = fullfile(outputDir, 'globalBeta_outputs');
    end

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
end


function outFile = output_file_in_dir(cfg, configuredFile, suffix, outputDir)
    outputFile = configured_or_input_output_file(cfg, configuredFile, suffix);
    outFile = fullfile_if_no_folder(outputDir, outputFile);
end


function outFile = fullfile_if_no_folder(outputDir, outputFile)
    outputFile = char(outputFile);
    [folderPart, ~, ~] = fileparts(outputFile);

    if isempty(folderPart)
        outFile = fullfile(outputDir, outputFile);
    else
        outFile = outputFile;
    end
end


function outFile = configured_or_input_output_file(cfg, configuredFile, suffix)
    if is_empty_text(configuredFile)
        outFile = input_based_output_file(cfg, suffix);
    else
        outFile = char(configuredFile);
    end
end


function outFile = input_based_output_file(cfg, suffix)
    outFile = sprintf('%s_%s.mat', input_base_name(cfg), suffix);
end


function prefix = input_based_output_prefix(cfg, suffix)
    prefix = sprintf('%s_%s', input_base_name(cfg), suffix);
end


function baseName = input_base_name(cfg)
    matFile = cfg.matFile;

    if isstring(matFile)
        matFile = char(matFile);
    elseif iscell(matFile)
        matFile = char(matFile{1});
    end

    [~, baseName, ~] = fileparts(matFile);

    if isempty(baseName)
        baseName = 'input';
    end

    baseName = regexprep(baseName, '[^a-zA-Z0-9_\-]+', '_');
    baseName = regexprep(baseName, '_+', '_');
    baseName = regexprep(baseName, '^_|_$', '');

    if isempty(baseName)
        baseName = 'input';
    end
end


function tf = is_empty_text(x)
    if isempty(x)
        tf = true;
    elseif isstring(x)
        tf = all(strlength(x) == 0);
    elseif ischar(x)
        tf = isempty(strtrim(x));
    else
        tf = false;
    end
end


function suffix = mode_index_suffix(modeIdx)
    modeIdx = modeIdx(:).';

    if isempty(modeIdx)
        suffix = 'modes';
        return;
    end

    suffix = ['mode' strjoin(compose('%d', modeIdx), '_')];
end


function pair = build_selected_mode_pair_reconstruction( ...
    q, t, hodmd, Uq, A_scale, cfg, r, d, Nt, modeTable)
    pair = struct();
    pair.valid = false;
    pair.A = [];
    pair.modeIdx = [];
    pair.representativeMode = [];
    pair.period_fs = nan;
    pair.decayTime_fs = nan;
    pair.wavenumber_cm1 = nan;
    pair.frequency_fs = nan;
    pair.growthRate_fs = nan;
    pair.qSpectrum = struct('valid', false);

    if ~isfield(cfg, 'modePairSave') || ...
            ~isfield(cfg.modePairSave, 'enabled') || ...
            ~cfg.modePairSave.enabled
        return;
    end

    modeIdx = cfg.modePairSave.modeIdx;
    modeIdx = unique_valid_indices(modeIdx, numel(hodmd.lambda));

    if isempty(modeIdx)
        return;
    end

    B_selected = reconstruct_hodmd_modes(hodmd, modeIdx, r, d, Nt);
    A_selectedModes = Uq * B_selected;
    A_selectedModes = undo_component_scale(A_selectedModes, A_scale, cfg);

    repMode = choose_representative_mode(modeIdx, modeTable);

    pair.valid = true;
    pair.A = A_selectedModes;
    pair.modeIdx = modeIdx;
    pair.representativeMode = repMode;
    pair.period_fs = modeTable.Period_fs(repMode);
    pair.decayTime_fs = modeTable.DecayTime_fs(repMode);
    pair.wavenumber_cm1 = modeTable.Wavenumber_cm1(repMode);
    pair.frequency_fs = modeTable.Frequency(repMode);
    pair.growthRate_fs = modeTable.GrowthRate(repMode);
    pair.qSpectrum = compute_hodmd_q_spectrum( ...
        q, t, A_selectedModes, hodmd, modeIdx);
end


function saveInfo = save_selected_mode_pair_reconstruction_if_requested( ...
    q, t, hodmd, Uq, A_scale, cfg, r, d, Nt, modeTable)
    saveInfo = struct();
    saveInfo.enabled = false;
    saveInfo.saved = false;
    saveInfo.file = '';
    saveInfo.modeIdx = [];
    saveInfo.representativeMode = [];
    saveInfo.period_fs = nan;
    saveInfo.decayTime_fs = nan;
    saveInfo.wavenumber_cm1 = nan;

    if ~isfield(cfg, 'modePairSave') || ...
            ~isfield(cfg.modePairSave, 'enabled') || ...
            ~cfg.modePairSave.enabled
        return;
    end

    saveInfo.enabled = true;

    pair = build_selected_mode_pair_reconstruction( ...
        q, t, hodmd, Uq, A_scale, cfg, r, d, Nt, modeTable);

    if ~pair.valid
        warning('Mode-pair save skipped because cfg.modePairSave.modeIdx is empty or invalid.');
        return;
    end

    A_selectedModes = pair.A;
    A_mode23 = A_selectedModes;
    A_mode23_qDeltaS = A_mode23;

    period_fs = pair.period_fs;
    decayTime_fs = pair.decayTime_fs;
    wavenumber_cm1 = pair.wavenumber_cm1;
    frequency_fs = pair.frequency_fs;
    growthRate_fs = pair.growthRate_fs;
    selectedModeIdx = pair.modeIdx;
    representativeMode = pair.representativeMode;
    qSpectrum = pair.qSpectrum;
    qSpectrumCombined = qSpectrum.qCombined;
    qSpectrumModes = qSpectrum.qModes;
    qSpectrumBasisCoefficients = qSpectrum.C;
    qSpectrumBasisInfo = qSpectrum.basisInfo;
    qSpectrumFitBasis = qSpectrum.G;
    qSpectrumMethod = qSpectrum.method;
    qSpectrumCondG = qSpectrum.condG;
    qSpectrumRankG = qSpectrum.rankG;
    qSpectrumModeReal = qSpectrum.qModeReal;
    qSpectrumModeImag = qSpectrum.qModeImag;
    qSpectrumModeAbs = qSpectrum.qModeAbs;
    qSpectrumCos = qSpectrum.C_cos;
    qSpectrumSin = qSpectrum.C_sin;
    qSpectrumAmp = qSpectrum.C_amp;
    qSpectrumPhase = qSpectrum.C_phase;

    outputDir = global_beta_output_dir(cfg);
    if isfield(cfg.modePairSave, 'outputDir') && ~isempty(cfg.modePairSave.outputDir)
        outputDir = cfg.modePairSave.outputDir;
    end

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    outputFile = input_based_output_file( ...
        cfg, ['globalBeta_' mode_index_suffix(cfg.modePairSave.modeIdx)]);
    if isfield(cfg.modePairSave, 'outputFile') && ~isempty(cfg.modePairSave.outputFile)
        outputFile = configured_or_input_output_file( ...
            cfg, cfg.modePairSave.outputFile, ...
            ['globalBeta_' mode_index_suffix(cfg.modePairSave.modeIdx)]);
    end

    outFile = fullfile_if_no_folder(outputDir, outputFile);

    save(outFile, ...
        'A_mode23', 'A_mode23_qDeltaS', 'A_selectedModes', ...
        'q', 't', ...
        'selectedModeIdx', 'representativeMode', ...
        'period_fs', 'decayTime_fs', ...
        'wavenumber_cm1', 'frequency_fs', 'growthRate_fs', ...
        'qSpectrum', 'qSpectrumCombined', 'qSpectrumModes', ...
        'qSpectrumBasisCoefficients', 'qSpectrumBasisInfo', ...
        'qSpectrumFitBasis', 'qSpectrumMethod', 'qSpectrumCondG', ...
        'qSpectrumRankG', ...
        'qSpectrumModeReal', 'qSpectrumModeImag', 'qSpectrumModeAbs', ...
        'qSpectrumCos', 'qSpectrumSin', 'qSpectrumAmp', 'qSpectrumPhase');
    append_input_fwhm_if_available(outFile, cfg);

    saveInfo.saved = true;
    saveInfo.file = outFile;
    saveInfo.modeIdx = selectedModeIdx;
    saveInfo.representativeMode = representativeMode;
    saveInfo.period_fs = period_fs;
    saveInfo.decayTime_fs = decayTime_fs;
    saveInfo.wavenumber_cm1 = wavenumber_cm1;
    saveInfo.frequency_fs = frequency_fs;
    saveInfo.growthRate_fs = growthRate_fs;
    saveInfo.qSpectrum = qSpectrum;

    fprintf('\nSaved selected mode-pair 2D qDeltaS reconstruction:\n');
    fprintf('  modes = %s\n', mat2str(selectedModeIdx));
    fprintf('  representative mode = %d\n', representativeMode);
    fprintf('  period = %.6g fs\n', period_fs);
    fprintf('  damping time constant = %.6g fs\n', decayTime_fs);
    fprintf('  file = %s\n\n', outFile);
end


function save_selected_mode_pair_batch_if_requested(selectedModePairBatch, cfg)
    if ~isfield(cfg, 'modePairSave') || ...
            ~isfield(cfg.modePairSave, 'enabled') || ...
            ~cfg.modePairSave.enabled
        return;
    end

    outputDir = global_beta_output_dir(cfg);
    if isfield(cfg.modePairSave, 'outputDir') && ~isempty(cfg.modePairSave.outputDir)
        outputDir = cfg.modePairSave.outputDir;
    end

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    outputFile = input_based_output_file( ...
        cfg, ['globalBeta_' mode_index_suffix(cfg.modePairSave.modeIdx)]);
    if isfield(cfg.modePairSave, 'outputFile') && ~isempty(cfg.modePairSave.outputFile)
        outputFile = configured_or_input_output_file( ...
            cfg, cfg.modePairSave.outputFile, ...
            ['globalBeta_' mode_index_suffix(cfg.modePairSave.modeIdx)]);
    end

    outFile = fullfile_if_no_folder(outputDir, outputFile);

    A_mode23_all = selectedModePairBatch.A_mode23_all;
    A_mode23_qDeltaS_all = A_mode23_all;
    A_selectedModes_all = selectedModePairBatch.A_selectedModes_all;
    q = selectedModePairBatch.q;
    t = selectedModePairBatch.t;
    selectedModeIdx = selectedModePairBatch.selectedModeIdx;
    representativeMode = selectedModePairBatch.representativeMode;
    period_fs = selectedModePairBatch.period_fs;
    decayTime_fs = selectedModePairBatch.decayTime_fs;
    wavenumber_cm1 = selectedModePairBatch.wavenumber_cm1;
    frequency_fs = selectedModePairBatch.frequency_fs;
    growthRate_fs = selectedModePairBatch.growthRate_fs;
    qSpectrumCombined_all = selectedModePairBatch.qSpectrumCombined_all;
    qSpectrumModes_all = selectedModePairBatch.qSpectrumModes_all;
    qSpectrumBasisCoefficients_all = selectedModePairBatch.qSpectrumBasisCoefficients_all;
    qSpectrumModeReal_all = selectedModePairBatch.qSpectrumModeReal_all;
    qSpectrumModeImag_all = selectedModePairBatch.qSpectrumModeImag_all;
    qSpectrumModeAbs_all = selectedModePairBatch.qSpectrumModeAbs_all;
    qSpectrumCos_all = selectedModePairBatch.qSpectrumCos_all;
    qSpectrumSin_all = selectedModePairBatch.qSpectrumSin_all;
    qSpectrumAmp_all = selectedModePairBatch.qSpectrumAmp_all;
    qSpectrumPhase_all = selectedModePairBatch.qSpectrumPhase_all;
    qSpectrumNumBasis = selectedModePairBatch.qSpectrumNumBasis;
    qSpectrumCondG = selectedModePairBatch.qSpectrumCondG;
    qSpectrumRankG = selectedModePairBatch.qSpectrumRankG;
    qSpectrumBasisInfo = selectedModePairBatch.qSpectrumBasisInfo;

    qSpectrumBatch = struct();
    qSpectrumBatch.q = q;
    qSpectrumBatch.selectedModeIdx = selectedModeIdx;
    qSpectrumBatch.qCombined_all = qSpectrumCombined_all;
    qSpectrumBatch.qModes_all = qSpectrumModes_all;
    qSpectrumBatch.qBasisCoefficients_all = qSpectrumBasisCoefficients_all;
    qSpectrumBatch.qModeReal_all = qSpectrumModeReal_all;
    qSpectrumBatch.qModeImag_all = qSpectrumModeImag_all;
    qSpectrumBatch.qModeAbs_all = qSpectrumModeAbs_all;
    qSpectrumBatch.qCos_all = qSpectrumCos_all;
    qSpectrumBatch.qSin_all = qSpectrumSin_all;
    qSpectrumBatch.qAmp_all = qSpectrumAmp_all;
    qSpectrumBatch.qPhase_all = qSpectrumPhase_all;
    qSpectrumBatch.numBasis = qSpectrumNumBasis;
    qSpectrumBatch.condG = qSpectrumCondG;
    qSpectrumBatch.rankG = qSpectrumRankG;
    qSpectrumBatch.basisInfo = qSpectrumBasisInfo;
    qSpectrumBatch.method = 'dehankelized_reconstruction_temporal_regression';

    save(outFile, ...
        'A_mode23_all', 'A_mode23_qDeltaS_all', 'A_selectedModes_all', ...
        'q', 't', ...
        'selectedModeIdx', 'representativeMode', ...
        'period_fs', 'decayTime_fs', ...
        'wavenumber_cm1', 'frequency_fs', 'growthRate_fs', ...
        'qSpectrumBatch', ...
        'qSpectrumCombined_all', 'qSpectrumModes_all', ...
        'qSpectrumBasisCoefficients_all', ...
        'qSpectrumModeReal_all', 'qSpectrumModeImag_all', ...
        'qSpectrumModeAbs_all', ...
        'qSpectrumCos_all', 'qSpectrumSin_all', ...
        'qSpectrumAmp_all', 'qSpectrumPhase_all', ...
        'qSpectrumNumBasis', 'qSpectrumCondG', 'qSpectrumRankG', ...
        'qSpectrumBasisInfo', ...
        '-v7.3');
    append_input_fwhm_if_available(outFile, cfg);

    fprintf('\nSaved selected mode-pair 3D qDeltaS reconstructions:\n');
    fprintf('  modes = %s\n', mat2str(cfg.modePairSave.modeIdx));
    fprintf('  file = %s\n\n', outFile);
end


function out = estimate_i2_distance_from_target_modes(t, q, A_target, modeIdx, cfg)
    % Estimate I-I distance from target-mode oscillatory q signal.
    %
    % The mock data generator uses R(t) around R_eq. Target HODMD modes
    % represent the oscillatory component, so this inverse uses the
    % linearized Debye response d(qDeltaS)/dR at R_eq.

    p = cfg.mockI2Distance;
    t = t(:);
    q = q(:);

    if size(A_target, 1) ~= numel(q) || size(A_target, 2) ~= numel(t)
        error('A_target must be Nq x Nt for I2 distance estimation.');
    end

    R_true = p.R_eq_A + p.amp_A .* exp(-t ./ p.tau_fs) .* ...
        cos(2*pi .* t ./ p.period_fs - p.phi_rad);

    fitBasis = i2_debye_dqdeltaS_dR(q, p.R_eq_A);

    if isfield(cfg.pre, 'smoothAlongQ') && cfg.pre.smoothAlongQ
        wq = round(cfg.pre.qSmoothWindow);
        wq = max(1, min(wq, numel(q)));
        fitBasis = movmean(fitBasis, wq, 1, 'Endpoints','shrink');
    end

    denom = fitBasis' * fitBasis;
    if denom <= eps
        error('I2 distance inverse failed because the Debye derivative basis is nearly zero.');
    end

    deltaR_est = (fitBasis' * A_target) ./ denom;
    deltaR_est = deltaR_est(:);
    R_est = p.R_eq_A + deltaR_est;

    A_fit = fitBasis * deltaR_est.';
    A_residual = A_target - A_fit;

    valid = isfinite(R_true) & isfinite(R_est);
    if any(valid)
        rmsError = sqrt(mean((R_est(valid) - R_true(valid)).^2));
    else
        rmsError = nan;
    end

    out = struct();
    out.t_fs = t;
    out.q = q;
    out.true_R_A = R_true;
    out.estimated_R_A = R_est;
    out.true_deltaR_A = R_true - p.R_eq_A;
    out.estimated_deltaR_A = deltaR_est;
    out.fitBasis_q = fitBasis;
    out.fittedTargetSignal = A_fit;
    out.fitResidual = A_residual;
    out.targetModeIdx = modeIdx;
    out.R_eq_A = p.R_eq_A;
    out.rmsError_A = rmsError;

    fprintf('I2 distance comparison from target modes:\n');
    fprintf('  target modes = %s\n', mat2str(modeIdx));
    fprintf('  linearized Debye basis around R_eq = %.4f Angstrom\n', p.R_eq_A);
    fprintf('  RMS error vs generator R(t) = %.6g Angstrom\n\n', rmsError);
end


function J = i2_debye_dqdeltaS_dR(q, R)
    q = q(:);
    fI = i2_atomic_form_factor(q);

    x = q .* R;
    dSinXOverX = (x .* cos(x) - sin(x)) ./ (x.^2);

    small = abs(x) < sqrt(eps);
    dSinXOverX(small) = 0;

    J = q .* (2 .* fI.^2 .* q .* dSinXOverX);
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

function plot_qt_map(t, q, A, ttl, cfg)
    figure('Color','w','Name',ttl, ...
        'WindowStyle','docked');

    if size(A, 1) ~= numel(q) || size(A, 2) ~= numel(t)
        error('plot_qt_map expects A to be Nq x Nt.');
    end

    imagesc(q, t, A.');
    set(gca, 'YDir', 'reverse');
    xlabel('q');
    ylabel('time [fs]');
    title(ttl, 'Interpreter','none');
    colorbar;
    colormap(turbo);
    axis square;

    if cfg.plot.useSymmetricColorLimit
        c = prctile(abs(A(:)), cfg.plot.percentileForCLim);
        if c > 0
            clim([-c c]);
        end
    end

    axis square;
end


function plot_all_mode_component_reconstructions(t, q, modeComponentRecon, cfg)
    components = modeComponentRecon.components;
    nComp = numel(components);

    if nComp == 0
        return;
    end

    nCols = min(3, nComp);
    nRows = ceil(nComp / nCols);

    figure('Color','w','Name','HODMD all mode 2D reconstructions', ...
        'WindowStyle','docked');

    tiledlayout(nRows, nCols, 'TileSpacing','compact', 'Padding','compact');

    for k = 1:nComp
        A = components(k).A;

        if size(A, 1) ~= numel(q) || size(A, 2) ~= numel(t)
            error('Component map size must be Nq x Nt.');
        end

        nexttile;
        imagesc(q, t, A.');
        set(gca, 'YDir','reverse');
        xlabel('q');
        ylabel('time [fs]');
        title(mode_component_title(components(k)), 'Interpreter','none');
        colorbar;
        colormap(turbo);

        % Scale every component independently so weak modes remain visible.
        clim(compute_qt_clim(A(:), cfg));

        axis square;
    end
end


function titleText = mode_component_title(component)
    if component.isComplex
        if strcmp(component.type, 'complex_pair')
            kindText = 'complex pair';
        else
            kindText = 'complex unpaired';
        end
    else
        kindText = 'real mode';
    end

    if isfinite(component.period_fs)
        periodText = sprintf('%.1f fs, %.1f cm^{-1}', ...
            component.period_fs, component.wavenumber_cm1);
    else
        periodText = sprintf('zero frequency, %.1f cm^{-1}', ...
            component.wavenumber_cm1);
    end

    if isfinite(component.decayTime_fs)
        if component.decayTime_fs > 0
            decayText = sprintf('decay time = %.1f fs', ...
                component.decayTime_fs);
        elseif component.decayTime_fs < 0
            decayText = sprintf('growth time = %.1f fs', ...
                abs(component.decayTime_fs));
        else
            decayText = 'decay time = Inf';
        end
    else
        decayText = 'decay time = Inf';
    end

    titleText = {sprintf('%s: %s, eigencomponents %s', ...
        component.dynamicMode, kindText, mat2str(component.modeIdx)), ...
        periodText, sprintf('%s, amplitude = %.5g', decayText, component.amplitude)};
end


function plot_i2_distance_comparison(distInfo)
    figure('Color','w','Name','I2 distance from target HODMD modes', ...
        'WindowStyle','docked');

    plot(distInfo.t_fs, distInfo.true_R_A, 'k-', ...
        'LineWidth', 1.8, 'DisplayName','generator R(t)');
    hold on;
    plot(distInfo.t_fs, distInfo.estimated_R_A, 'r--', ...
        'LineWidth', 1.6, 'DisplayName','target-mode inverse R(t)');
    yline(distInfo.R_eq_A, 'k:', ...
        sprintf('R_{eq}=%.3f', distInfo.R_eq_A), ...
        'DisplayName','R_eq');
    hold off;

    xlabel('time [fs]');
    ylabel('I-I distance [Angstrom]');
    title(sprintf('I2 distance from target modes %s, RMS error %.4g Angstrom', ...
        mat2str(distInfo.targetModeIdx), distInfo.rmsError_A), ...
        'Interpreter','none');
    legend('Location','best');
    grid on;
end


function plot_batch_target_statistics(batchResults)
    T = batchResults.summaryTable;
    x = batchResults.targetWavenumber_cm1;

    if ~isfield(batchResults, 'truth')
        batchResults.truth = struct('enabled', false, ...
            'period_fs', nan, 'wavenumber_cm1', nan, 'dampingTime_fs', nan);
    end

    figure('Color','w','Name','3D batch nearest-wavenumber HODMD statistics', ...
        'WindowStyle','docked');

    tiledlayout(1, 3, 'TileSpacing','compact', 'Padding','compact');

    nexttile;
    plot_batch_raw_points(x, batchResults.Period_fs);
    hold on;
    plot_batch_truth_point(x, batchResults.truth, batchResults.truth.period_fs);
    errorbar(x, T.PeriodMean_fs, T.PeriodStd_fs, ...
        'ko-', 'LineWidth', 1.3, 'MarkerFaceColor','k');
    hold off;
    ylabel('period [fs]');
    title('Selected-mode period: individual runs and mean +/- std');
    grid on;
    set_batch_target_xlim(x);
    add_batch_stats_text(x, T.PeriodMean_fs, T.PeriodStd_fs, 'fs', ...
        sum(isfinite(batchResults.Period_fs), 1), batchResults.nSeries);
    hide_batch_x_axis();

    nexttile;
    plot_batch_raw_points(x, batchResults.Wavenumber_cm1);
    hold on;
    plot_batch_truth_point(x, batchResults.truth, batchResults.truth.wavenumber_cm1);
    errorbar(x, T.WavenumberMean_cm1, T.WavenumberStd_cm1, ...
        'ko-', 'LineWidth', 1.3, 'MarkerFaceColor','k');
    hold off;
    ylabel('wavenumber [cm^{-1}]');
    title('Selected-mode wavenumber: individual runs and mean +/- std');
    grid on;
    set_batch_target_xlim(x);
    add_batch_stats_text(x, T.WavenumberMean_cm1, T.WavenumberStd_cm1, ...
        'cm^{-1}', sum(isfinite(batchResults.Wavenumber_cm1), 1), batchResults.nSeries);
    hide_batch_x_axis();

    nexttile;
    plot_batch_raw_points(x, batchResults.DecayTime_fs);
    hold on;
    plot_batch_truth_point(x, batchResults.truth, batchResults.truth.dampingTime_fs);
    errorbar(x, T.DecayTimeMean_fs, T.DecayTimeStd_fs, ...
        'ko-', 'LineWidth', 1.3, 'MarkerFaceColor','k');
    hold off;
    ylabel('damping time constant [fs]');
    title('Selected-mode damping time constant: individual runs and mean +/- std');
    grid on;
    set_batch_target_xlim(x);
    add_batch_stats_text(x, T.DecayTimeMean_fs, T.DecayTimeStd_fs, 'fs', ...
        sum(isfinite(batchResults.DecayTime_fs), 1), batchResults.nSeries);
    hide_batch_x_axis();
    sgtitle(sprintf('Nearest pair to %.4g cm^{-1}: %d / %d runs (no tolerance limit)', ...
        x, T.NSelected, batchResults.nSeries));
end


function add_batch_stats_text(targetPeriods, meanVals, stdVals, unitText, ...
        nValid, nSeries)
    targetPeriods = targetPeriods(:);
    meanVals = meanVals(:);
    stdVals = stdVals(:);
    nValid = nValid(:);

    valid = isfinite(targetPeriods) & isfinite(meanVals) & ...
        isfinite(stdVals) & isfinite(nValid);

    if ~any(valid)
        return;
    end

    targetPeriods = targetPeriods(valid);
    meanVals = meanVals(valid);
    stdVals = stdVals(valid);
    nValid = nValid(valid);

    nTargets = numel(targetPeriods);
    lines = strings(2*nTargets, 1);

    for iTarget = 1:nTargets
        lines(2*iTarget-1) = compose( ...
            'target %.4g cm^{-1}: %.4g +/- %.4g %s', ...
            targetPeriods(iTarget), meanVals(iTarget), ...
            stdVals(iTarget), unitText);
        lines(2*iTarget) = compose( ...
            'finite statistics: %d / %d', ...
            round(nValid(iTarget)), nSeries);
    end

    text(0.02, 0.96, lines, ...
        'Units','normalized', ...
        'VerticalAlignment','top', ...
        'HorizontalAlignment','left', ...
        'BackgroundColor','w', ...
        'EdgeColor', [0.82 0.82 0.82], ...
        'Margin', 4, ...
        'Interpreter','tex');
end


function hide_batch_x_axis()
    ax = gca;
    ax.XTick = [];
    ax.XColor = 'none';
    ax.XLabel.String = '';
end


function plot_batch_truth_point(x, truth, yTruth)
    if ~isfield(truth, 'enabled') || ~truth.enabled || ...
            ~isfinite(truth.period_fs) || ~isfinite(yTruth)
        return;
    end

    x = x(:);

    if isempty(x)
        return;
    end

    [~, idx] = min(abs(x - truth.period_fs));
    xTruth = x(idx);

    scatter(xTruth, yTruth, 70, ...
        'MarkerFaceColor', 'r', ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 0.8);
end


function plot_batch_raw_points(x, Y)
    x = x(:).';
    nTargets = numel(x);

    if isempty(Y)
        return;
    end

    hold on;

    for iTarget = 1:nTargets
        vals = Y(:, iTarget);
        vals = vals(isfinite(vals));

        if isempty(vals)
            continue;
        end

        xj = x(iTarget) .* ones(size(vals));

        scatter(xj, vals, 12, ...
            'MarkerFaceColor', [0.25 0.45 0.90], ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.35);
    end
end


function set_batch_target_xlim(x)
    x = x(isfinite(x));

    if isempty(x)
        return;
    end

    if isscalar(x) || max(x) == min(x)
        pad = max(10, 0.1 * abs(x(1)));
        xlim([x(1)-pad, x(1)+pad]);
    else
        pad = 0.08 * (max(x) - min(x));
        xlim([min(x)-pad, max(x)+pad]);
    end
end


function climVals = compute_qt_clim(vals, cfg)
    vals = vals(isfinite(vals));

    if isempty(vals)
        climVals = [-1 1];
    elseif cfg.plot.useSymmetricColorLimit
        c = prctile(abs(vals), cfg.plot.percentileForCLim);
        if ~isfinite(c) || c <= 0
            c = max(abs(vals));
        end
        if ~isfinite(c) || c <= 0
            c = 1;
        end
        climVals = [-c c];
    else
        lo = min(vals);
        hi = max(vals);
        if lo == hi
            pad = max(abs(lo), 1) * 0.01;
            lo = lo - pad;
            hi = hi + pad;
        end
        climVals = [lo hi];
    end
end


function markerSize = log_scaled_marker_size(values, minSize, maxSize)
    values = abs(values(:));
    markerSize = minSize * ones(size(values));

    finitePositive = isfinite(values) & values > 0;

    if ~any(finitePositive)
        return;
    end

    minPositive = min(values(finitePositive));
    safeValues = values;
    safeValues(~finitePositive) = minPositive;

    logValues = log10(safeValues);
    logMin = min(logValues);
    logMax = max(logValues);

    if logMax <= logMin
        markerSize(:) = 0.5 * (minSize + maxSize);
        return;
    end

    markerSize = minSize + ...
        (maxSize - minSize) .* (logValues - logMin) ./ (logMax - logMin);
end


function plot_hodmd_diagnostics(hodmd, dmdInfo, modeTable, cfg)
    % 1. Delay SVD spectrum
    s = dmdInfo.delaySVDValues(:);
    ce = dmdInfo.delayCumulativeEnergy(:);
    sFrac = s / max(sum(s), eps);

    rShow = min(cfg.plot.maxModesToShow, numel(s));

    if isfield(dmdInfo, 'isDelayEmbedded') && ~dmdInfo.isDelayEmbedded
        figName = 'DMD input SVD spectrum';
        spectrumTitle = 'SVD spectrum of direct DMD snapshots Z1';
        cumulativeTitle = 'DMD input SVD cumulative energy';
    else
        figName = 'Delay-embedded SVD spectrum';
        spectrumTitle = 'SVD spectrum of delay-embedded snapshots';
        cumulativeTitle = 'Delay SVD cumulative energy';
    end

    figure('Color','w','Name',figName);

    subplot(2,1,1);
    semilogy(1:rShow, sFrac(1:rShow), 'o-', 'LineWidth', 1.2);
    xlabel('delay SVD index');
    ylabel('singular value fraction');
    title(spectrumTitle);
    grid on;

    subplot(2,1,2);
    plot(1:rShow, ce(1:rShow), 'o-', 'LineWidth', 1.2);
    yline(0.90, '--');
    yline(0.95, '--');
    yline(0.99, '--');
    xlabel('delay SVD index');
    ylabel('cumulative energy');
    title(cumulativeTitle);
    ylim([0 1.02]);
    grid on;


    % 2. DMD eigenvalues
    lambda = hodmd.lambda(:);
    amp = modeTable.Amplitude;
    markerSizeEigen = log_scaled_marker_size(amp, 20, 140);
    markerSize = 20 + 80 * amp / max(amp + eps);

    figure('Color','w','Name','HODMD eigenvalues');

    th = linspace(0, 2*pi, 400);
    plot(cos(th), sin(th), 'k--', 'LineWidth', 1.0);
    hold on;
    scatter(real(lambda), imag(lambda), markerSizeEigen, 'filled');
    axis equal;
    grid on;
    xlabel('real(\lambda)');
    ylabel('imag(\lambda)');
    title('HODMD eigenvalues in complex plane, marker size = log_{10}(amplitude)');
    hold off;


    % 3. Period vs damping time constant
    finitePeriod = isfinite(modeTable.Period_fs) & modeTable.Frequency > 0;
    decayMode = strcmp(modeTable.Stability, 'decay') & ...
        isfinite(modeTable.DecayTime_fs) & modeTable.DecayTime_fs > 0;
    plotMask = finitePeriod & decayMode;

    figure('Color','w','Name','HODMD period vs damping time');

    scatter(modeTable.Period_fs(plotMask), ...
            modeTable.DecayTime_fs(plotMask), ...
            markerSize(plotMask), 'filled');

    hold on;
    for i = 1:numel(cfg.target.period_fs)
        xline(cfg.target.period_fs(i), '--', sprintf('%.0f fs', cfg.target.period_fs(i)));
    end
    hold off;

    xlabel('period [fs]');
    ylabel('damping time constant [fs]');
    title('HODMD modes: period vs damping time constant');
    grid on;

    xlim([0 max([cfg.target.period_fs(:)'*1.8, 500])]);
end


function inspect_reconstructed_signal(Arec, t, cfg, labelText)
    % Representative temporal score by SVD of Arec(q,t),
    % plus FFT of that score.

    [~, S, V] = svd(Arec, 'econ');

    score = S(1,1) * V(:,1);
    score = score(:).';

    [~, imax] = max(abs(score));
    if score(imax) < 0
        score = -score;
    end

    dt = median(diff(t));

    [freq, pow] = onesided_fft_power(score - mean(score), dt);

    figure('Color','w','Name',['Temporal inspection - ' labelText]);

    subplot(2,1,1);
    plot(t, score, 'LineWidth', 1.2);
    xlabel('time [fs]');
    ylabel('representative temporal score');
    title([labelText ' temporal score'], 'Interpreter','none');
    grid on;

    subplot(2,1,2);
    plot(freq, pow / max(pow + eps), 'LineWidth', 1.2);
    xlabel('frequency [cycles/fs]');
    ylabel('normalized FFT power');
    title([labelText ' temporal FFT'], 'Interpreter','none');
    grid on;

    if isfield(cfg.plot, 'maxFreqToPlot') && ~isempty(cfg.plot.maxFreqToPlot)
        xlim([0 cfg.plot.maxFreqToPlot]);
    end

    hold on;
    for i = 1:numel(cfg.target.period_fs)
        f0 = 1 / cfg.target.period_fs(i);
        xline(f0, '--', sprintf('%.0f fs', cfg.target.period_fs(i)));
    end
    hold off;
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


function qSpectrum = plot_hodmd_q_spectrum(q, t, Uq, hodmd, modeIdx, r, d, A_scale, cfg, labelText)
    % Plot q-spectrum of selected HODMD/DMD mode(s).
    %
    % The selected mode trajectory is reconstructed in delay space,
    % de-Hankelized back to the original time series, lifted to q-space,
    % then each q trace is fit to the modal temporal basis.

    B_rec = reconstruct_hodmd_modes(hodmd, modeIdx, r, d, numel(t));
    A_rec = Uq * B_rec;
    A_rec = undo_component_scale(A_rec, A_scale, cfg);

    qSpectrum = compute_hodmd_q_spectrum(q, t, A_rec, hodmd, modeIdx);

    if ~qSpectrum.valid
        return;
    end

    modeIdx = qSpectrum.modeIdx;
    qModes = qSpectrum.C;
    qCombined = qSpectrum.qCombined;
    basisLabels = q_spectrum_basis_labels(qSpectrum.basisInfo);

    save_hodmd_q_spectrum_if_requested(qSpectrum, modeIdx, cfg, labelText);

    figure('Color','w','Name',['HODMD q-spectrum - ' labelText]);

    subplot(2,1,1);
    plot(q, qCombined, 'LineWidth', 1.4);
    yline(0, 'k--');
    xlabel('q');
    ylabel('combined q spectrum');
    title([labelText ' combined q-spectrum'], 'Interpreter','none');
    grid on;

    subplot(2,1,2);
    hold on;
    for j = 1:size(qModes, 2)
        plot(q, qModes(:, j), 'LineWidth', 1.1);
    end
    hold off;
    xlabel('q');
    ylabel('basis coefficient');
    title([labelText ' modal-basis regression q coefficients'], 'Interpreter','none');
    legend(basisLabels, 'Location','best');
    grid on;
end


function spectrum = plot_hodmd_r_spectrum(qSpectrum, settings, labelText)
    spectrum = struct('valid', false);
    if ~qSpectrum.valid
        return;
    end
    validateattributes(settings.range_A, {'numeric'}, ...
        {'real', 'finite', 'nonnegative', 'vector', 'numel', 2});
    validateattributes(settings.step_A, {'numeric'}, ...
        {'real', 'finite', 'positive', 'scalar'});
    assert(settings.range_A(2) > settings.range_A(1), ...
        'cfg.rSpectrum.range_A must be an increasing range.');
    validateattributes(settings.inputIsQWeighted, {'logical', 'numeric'}, ...
        {'scalar', 'binary'});

    [q, order] = sort(qSpectrum.q(:));
    assert(numel(q) >= 2 && all(isfinite(q)) && all(diff(q) > 0), ...
        'Sine transform requires at least two distinct finite q points.');
    C = qSpectrum.C(order, :);
    integrandCoefficients = C;
    if ~settings.inputIsQWeighted
        integrandCoefficients = q .* C;
    end
    r = (settings.range_A(1):settings.step_A:settings.range_A(2)).';
    assert(numel(r) >= 2, 'The r range must contain at least two points.');

    % Trapezoidal quadrature also supports nonuniform q sampling.
    dq = diff(q);
    weights = [dq(1); dq(1:end-1) + dq(2:end); dq(end)] / 2;
    F = (2/pi) * sin(r * q.') * (weights .* integrandCoefficients);

    % At the temporal-basis origin, cos/real weights are 1 and sin weights 0.
    if isempty(settings.evaluationTime_fs)
        settings.evaluationTime_fs = qSpectrum.t(1);
    end
    validateattributes(settings.evaluationTime_fs, {'numeric'}, ...
        {'real', 'finite', 'scalar'});
    % Match the unnormalized temporal basis used to fit qSpectrum.C.
    tau = settings.evaluationTime_fs - qSpectrum.t(1);
    temporalWeights = zeros(size(C, 2), 1);
    for j = 1:numel(temporalWeights)
        info = qSpectrum.basisInfo(j);
        envelope = exp(real(info.omega) * tau);
        phase = abs(imag(info.omega)) * tau;
        switch info.type
            case 'cos'
                temporalWeights(j) = envelope * cos(phase);
            case 'sin'
                temporalWeights(j) = envelope * sin(phase);
            case 'real'
                temporalWeights(j) = envelope;
            otherwise
                error('Unsupported modal basis type: %s', info.type);
        end
    end
    qAtTime = C * temporalWeights;
    rAtTime = F * temporalWeights;
    spectrum.valid = true;
    spectrum.r_A = r;
    spectrum.C = F;
    spectrum.q = q;
    spectrum.qCoefficients = C;
    spectrum.evaluationTime_fs = settings.evaluationTime_fs;
    spectrum.referenceTime_fs = qSpectrum.t(1);
    spectrum.temporalWeights = temporalWeights;
    spectrum.qAtTime = qAtTime;
    spectrum.rAtTime = rAtTime;
    spectrum.basisInfo = qSpectrum.basisInfo;
    spectrum.modeIdx = qSpectrum.modeIdx;
    spectrum.settings = settings;
    spectrum.method = 'finite_q_range_sine_transform_trapezoidal';

    figure('Color', 'w', 'WindowStyle', 'docked', ...
        'Name', ['HODMD target r-spectrum - ' labelText], ...
        'Tag', 'HODMDTargetRSpectrum');
    tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile;
    plot(q, qAtTime, 'LineWidth', 1.3);
    yline(0, 'k:');
    xlabel('q [Angstrom^{-1}]');
    ylabel('Target signal [a.u.]');
    title(sprintf('Target modes %s: q spectrum at t = %.6g fs', ...
        mat2str(qSpectrum.modeIdx), settings.evaluationTime_fs));
    grid on;
    nexttile;
    plot(r, rAtTime, 'LineWidth', 1.3);
    yline(0, 'k:');
    xlabel('r [Angstrom]');
    ylabel('Sine transform [a.u.]');
    title(sprintf('Sine transform at t = %.6g fs: q = %.3g to %.3g Angstrom^{-1}', ...
        settings.evaluationTime_fs, q(1), q(end)));
    xlim(settings.range_A);
    grid on;
end


function labels = q_spectrum_basis_labels(basisInfo)
    labels = cell(1, numel(basisInfo));

    for j = 1:numel(basisInfo)
        labels{j} = sprintf('%s mode %d', basisInfo(j).type, basisInfo(j).sourceMode);
    end
end


function save_hodmd_q_spectrum_if_requested(qSpectrum, modeIdx, cfg, labelText)
    % Save combined real q-spectrum and q to .mat if requested.
    %
    % Saved variables:
    %   q                  : q-axis, Nq x 1
    %   qCombined          : combined real q-spectrum, Nq x 1
    %   modeIdx            : selected HODMD/DMD mode indices
    %   labelText          : source reconstruction label
    %   qModes             : modal temporal-basis q coefficients

    if ~isfield(cfg, 'qSpectrum') || ...
            ~isfield(cfg.qSpectrum, 'saveMat') || ...
            ~cfg.qSpectrum.saveMat
        return;
    end

    outputDir = global_beta_output_dir(cfg);
    if isfield(cfg.qSpectrum, 'outputDir') && ~isempty(cfg.qSpectrum.outputDir)
        outputDir = cfg.qSpectrum.outputDir;
    end

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    filePrefix = input_based_output_prefix(cfg, 'globalBeta_qSpec');
    if isfield(cfg.qSpectrum, 'filePrefix') && ~isempty(cfg.qSpectrum.filePrefix)
        if is_empty_text(cfg.qSpectrum.filePrefix)
            filePrefix = input_based_output_prefix(cfg, 'globalBeta_qSpec');
        else
            filePrefix = char(cfg.qSpectrum.filePrefix);
        end
    end

    modeStr = strrep(mat2str(modeIdx), ' ', '_');
    modeStr = regexprep(modeStr, '[\[\]]', '');

    labelSafe = regexprep(labelText, '[^a-zA-Z0-9_\-]+', '_');
    labelSafe = regexprep(labelSafe, '_+', '_');
    labelSafe = regexprep(labelSafe, '^_|_$', '');

    if numel(labelSafe) > 80
        labelSafe = labelSafe(1:80);
    end

    outFile = fullfile(outputDir, sprintf('%s_%s_modes_%s.mat', ...
        filePrefix, labelSafe, modeStr));

    q = qSpectrum.q(:);
    qCombined = qSpectrum.qCombined(:);
    qModes = qSpectrum.qModes;
    qSpectrumBasisCoefficients = qSpectrum.C;
    qSpectrumBasisInfo = qSpectrum.basisInfo;
    qSpectrumFitBasis = qSpectrum.G;
    qSpectrumMethod = qSpectrum.method;
    qSpectrumCondG = qSpectrum.condG;
    qSpectrumRankG = qSpectrum.rankG;
    qSpectrumCos = qSpectrum.C_cos;
    qSpectrumSin = qSpectrum.C_sin;
    qSpectrumAmp = qSpectrum.C_amp;
    qSpectrumPhase = qSpectrum.C_phase;
    selectedModes = modeIdx(:).';

    saveIndividualModes = false;
    if isfield(cfg.qSpectrum, 'saveIndividualModes')
        saveIndividualModes = cfg.qSpectrum.saveIndividualModes;
    end

    if saveIndividualModes
        save(outFile, 'q', 'qCombined', 'qModes', ...
            'qSpectrumBasisCoefficients', 'qSpectrumBasisInfo', ...
            'qSpectrumFitBasis', 'qSpectrumMethod', 'qSpectrumCondG', ...
            'qSpectrumRankG', ...
            'qSpectrumCos', 'qSpectrumSin', 'qSpectrumAmp', 'qSpectrumPhase', ...
            'modeIdx', 'selectedModes', 'labelText');
    else
        save(outFile, 'q', 'qCombined', ...
            'qSpectrumBasisInfo', 'qSpectrumMethod', 'qSpectrumCondG', ...
            'qSpectrumRankG', ...
            'qSpectrumCos', 'qSpectrumSin', 'qSpectrumAmp', 'qSpectrumPhase', ...
            'modeIdx', 'selectedModes', 'labelText');
    end
    append_input_fwhm_if_available(outFile, cfg);

    fprintf('\nSaved combined q-spectrum to %s\n', outFile);
end


function append_input_fwhm_if_available(outFile, cfg)
    % Preserve input-file IRF metadata in every MAT result without changing units.

    if ~isfield(cfg, 'inputMetadata') || ...
            ~isfield(cfg.inputMetadata, 'hasFWHM') || ...
            ~cfg.inputMetadata.hasFWHM
        return;
    end

    FWHM = cfg.inputMetadata.FWHM;
    save(outFile, 'FWHM', '-append');
end


function [freq, pow] = onesided_fft_power(x, dt)
    x = x(:).';
    N = numel(x);

    x = x - mean(x);

    if N > 1
        w = 0.5 - 0.5*cos(2*pi*(0:N-1)/(N-1));  % Hann window
    else
        w = 1;
    end

    X = fft(x .* w);
    P2 = abs(X / N).^2;

    nKeep = floor(N/2) + 1;
    pow = P2(1:nKeep);

    if nKeep > 2
        pow(2:end-1) = 2 * pow(2:end-1);
    end

    freq = (0:nKeep-1) / (N * dt);
end


