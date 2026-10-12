% Step 01 (optional): generate I2 TRXL mock qDeltaS from the Debye equation.
% Outputs are written to data/mock in this distribution.

clear;
clc;

outDir = fullfile(fileparts(mfilename('fullpath')), 'data', 'mock');
if ~isfolder(outDir), mkdir(outDir); end

q = (1.5:0.02:6.5)';
t_fs = (300:20:2000)';
t = t_fs * 1e-15; % Save time delays in seconds in the MAT files.
irfDtFs = 1;
t_calc_fs = (0:irfDtFs:2300)';

noisePcts = 90;
rngSeed = 3; % Reproduces the archived noise90_seedN1 input; N1 is the count.
numRandomSeeds = 1;
irfFwhmFs = 100;

% I atomic form factor coefficients:
% f(q) = sum_i a_i * exp(-b_i * (q / 4pi)^2) + c
coeff = [20.1472 4.347 18.9949 0.3814 7.513 27.766 2.2735 66.8776 4.0712];
a = coeff([1 3 5 7]);
b = coeff([2 4 6 8]);
c = coeff(9);

x = (q ./ (4*pi)).^2;
fI = sum(a .* exp(-x * b), 2) + c;

% Table S1 damped cosine parameters for excited-state I-I distance.
R_eq = 3.300;
amp = 0.100;
wavenumber_cm = 110;
speedOfLight_cm_per_s = 2.99792458e10;
T = 1e15 / (speedOfLight_cm_per_s * wavenumber_cm);
tau = 400;
phi = pi;

% Ground-state I2 distance for difference scattering.
R_ground = 2.700;

R_excited = amp .* exp(-t_calc_fs' ./ tau) .* cos(2*pi .* t_calc_fs' ./ T - phi) + R_eq;

qR_excited = q * R_excited;
qR_ground = q .* R_ground;

S_excited = 2*fI.^2 + 2*fI.^2 .* sin(qR_excited) ./ qR_excited;
S_ground = 2*fI.^2 + 2*fI.^2 .* sin(qR_ground) ./ qR_ground;

qDeltaSRaw = q .* (S_excited - S_ground);

% Apply a normalized Gaussian IRF convolution along the time axis.
dtFs = mean(diff(t_calc_fs));
irfSigmaFs = irfFwhmFs / (2 * sqrt(2 * log(2)));
kernelHalfWidth = ceil(5 * irfSigmaFs / dtFs);
irfTimeFs = (-kernelHalfWidth:kernelHalfWidth) * dtFs;
irfKernel = exp(-0.5 * (irfTimeFs ./ irfSigmaFs).^2);
irfKernel = irfKernel ./ sum(irfKernel);

qDeltaSConv = conv2(qDeltaSRaw, irfKernel, 'same');
[isOutputTime, outputTimeIdx] = ismember(t_fs, t_calc_fs);
if ~all(isOutputTime)
    error('Output time grid must be a subset of the convolution time grid.');
end
qDeltaS = qDeltaSConv(:, outputTimeIdx);

% Noise definition: use the first right singular component of IRF-convolved qDeltaS.
[~, S, V] = svd(qDeltaS, 'econ');
rSV = S * V';
rSV1 = rSV(1, :);
baseSigma = (max(rSV1) - min(rSV1)) / 4;

baseNoise = zeros([size(qDeltaS), numRandomSeeds]);
for seedIdx = 1:numRandomSeeds
    rng(rngSeed + seedIdx - 1, 'twister');
    baseNoise(:, :, seedIdx) = randn(size(qDeltaS));
end

commonAbsMax = max(abs(qDeltaS(:)));
for k = 1:numel(noisePcts)
    noiseSigma = baseSigma * (noisePcts(k) / 100);
    thisAbsMax = max(abs(qDeltaS + noiseSigma .* baseNoise), [], 'all');
    commonAbsMax = max(commonAbsMax, thisAbsMax);
end

for k = 1:numel(noisePcts)
    noisePct = noisePcts(k);
    noiseSigma = baseSigma * (noisePct / 100);
    A = qDeltaS + noiseSigma .* baseNoise; %#ok<NASGU>

    filePrefix = sprintf('I2_qDeltaS_noise%d_seedN%d', noisePct, numRandomSeeds);
    matPath = fullfile(outDir, sprintf('%s.mat', filePrefix));
    save(matPath, 'A', 'q', 't');

    pngPath = fullfile(outDir, sprintf('%s_TRXL_square.png', filePrefix));

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 900]);
    ax = axes(fig);
    imagesc(ax, q, t_fs, A(:, :, 1)');
    set(ax, ...
        'YDir', 'reverse', ...
        'PlotBoxAspectRatio', [1 1 1], ...
        'CLim', [-commonAbsMax commonAbsMax]);
    colormap(ax, 'turbo');
    cb = colorbar(ax);
    cb.Label.String = 'q\DeltaS';
    xlabel(ax, 'q (A^{-1})');
    ylabel(ax, 'Time delay (fs)');
    title(ax, sprintf('I2 Debye q\\DeltaS, IRF %.0f fs FWHM, %.0f%% noise, seed 1/%d', ...
        irfFwhmFs, noisePct, numRandomSeeds));
    axis(ax, 'square');
    exportgraphics(fig, pngPath, 'Resolution', 180);
    close(fig);

    fprintf('noise %d%%: sigma = %.9g\n', noisePct, baseSigma * (noisePct / 100));
    fprintf('  saved %s\n', matPath);
    fprintf('  saved %s\n', pngPath);
end

fprintf('A size = %dx%d, q size = %dx%d, t size = %dx%d\n', ...
    size(qDeltaS, 1), size(qDeltaS, 2), size(q, 1), size(q, 2), size(t, 1), size(t, 2));
fprintf('Saved A size per MAT file = %dx%dx%d\n', size(qDeltaS, 1), size(qDeltaS, 2), numRandomSeeds);
fprintf('Random seed range = %d to %d\n', rngSeed, rngSeed + numRandomSeeds - 1);
fprintf('Saved t range = %.9g to %.9g seconds\n', t(1), t(end));
fprintf('IRF FWHM = %.9g fs, sigma = %.9g fs\n', irfFwhmFs, irfSigmaFs);
fprintf('Common color scale = +/- %.9g\n', commonAbsMax);
