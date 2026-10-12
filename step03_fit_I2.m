%% I2 90% noise, seedN1: direct oscillatory Debye fit only.
% Rgs/Req, IRF and absolute scale are fixed.
% Fit amplitude and phase only; T/tau are fixed to the selected HODMD eigenvalue.
% Only result files under this project's results directory are written.
clearvars; clc; close all;
set(groot,'defaultFigureWindowStyle','docked');
projectDir=fileparts(mfilename('fullpath'));
cfg.dataDir=fullfile(projectDir,'data','mock');
cfg.outputRoot=fullfile(projectDir,'results','I2_direct_noise90');
cfg.generatorFile=fullfile(projectDir,'step01_generate_noise.m');
cfg.sourceFile=fullfile(projectDir,'step02_run_HODMD.m');
cfg.noiseLevels=90;
cfg.seeds=1; % Fit the single realization in each seedN1 file.
cfg.fileSeedCount=1;
cfg.nStarts=3;
cfg.randomSeed=20261001;
cfg.Rgs=2.7; cfg.Req=3.3;
cfg.irfFwhm_fs=100; cfg.irfDt_fs=1; cfg.signalScale=1;
cfg.lower=[0.001 0]; cfg.upper=[0.3 2*pi];
cfg.initial=[0.08 2.5];
cfg.freeParameterIndices=[1 2];
cfg.fixedTimeSource='Selected positive-frequency HODMD eigenvalue';
cfg.maxIterations=100; cfg.maxFunctionEvaluations=1300;
cfg.functionTolerance=1e-8;cfg.stepTolerance=1e-7;cfg.optimalityTolerance=1e-6;
cfg.finiteDifferenceStep=1e-4;
cfg.makeFigures=true; cfg.exportFigures=true; cfg.runLabel='direct90_N1_fixed_T_tau';
cfg.hodmd=hodmd_core('defaults');
cfg.hodmd.pre.smoothAlongQ=true;cfg.hodmd.pre.qSmoothWindow=20;
% Raw I2 files are UNSMOOTHED. Apply MA20 once to each theory and raw input.
% Existing seedN100 results were numerically reproduced with this setting.
assert(cfg.Rgs==2.7 && cfg.Req==3.3,'This benchmark fixes Rgs=2.7 and Req=3.3.');
assert(exist('lsqnonlin','file')==2,'Optimization Toolbox is required.');
audit=hodmd_core('audit',cfg.sourceFile);assert(all(audit.MatchesOriginal));
cfg.generatorSource=fileread(cfg.generatorFile);
reference=[0.1 pi 1e15/(2.99792458e10*110) 400];
parameterNames={'Amplitude_A','Phase_rad','Period_fs','Decay_fs'};
rng(cfg.randomSeed,'twister');
starts=cfg.initial;
% Independent amplitude/phase starting values; no T/tau optimization coordinates.
for k=2:cfg.nStarts
    starts(k,:)=[0.04+0.12*rand,pi-1.2+2.4*rand];
end
assert(all(starts>=cfg.lower & starts<=cfg.upper,'all'),'Initial values exceed bounds.');
raw=cell(numel(cfg.noiseLevels),1);paths=raw;
for n=1:numel(cfg.noiseLevels)
    paths{n}=fullfile(cfg.dataDir,sprintf('I2_qDeltaS_noise%d_seedN%d.mat',cfg.noiseLevels(n),cfg.fileSeedCount));
    raw{n}=load(paths{n});assert(all(isfield(raw{n},{'A','q','t'})),'Missing A/q/t.');
    assert(ndims(raw{n}.A)<=3 && size(raw{n}.A,3)>=max(cfg.seeds),'Invalid seed selection.');
end
q=raw{1}.q(:); t=raw{1}.t(:).'*1e15;
for n=1:numel(raw)
    assert(isequal(raw{n}.q(:),q) && max(abs(raw{n}.t(:).'*1e15-t))<1e-8,'Data grids differ.');
    assert(size(raw{n}.A,1)==numel(q) && size(raw{n}.A,2)==numel(t),'Input orientation mismatch.');
    assert(all(isfinite(raw{n}.A(:))),'Nonfinite input data.');
end
c=struct('q',q,'t',t,'cfg',cfg,'f2',hodmd_core('formFactor',q).^2);
% Match the inspected generator EXACTLY: 0:1:2300 fs and a five-sigma kernel.
c.tf=0:cfg.irfDt_fs:2300;
sigma=cfg.irfFwhm_fs/(2*sqrt(2*log(2)));half=ceil(5*sigma/cfg.irfDt_fs);
k=exp(-0.5*((-half:half)*cfg.irfDt_fs/sigma).^2);c.kernel=k/sum(k);
c.sampleIdx=round(t/cfg.irfDt_fs)+1;
assert(max(abs(c.tf(c.sampleIdx)-t))<1e-8,'Output grid not on IRF grid.');
c.staticFine=2*q.*c.f2.*(sinc0(q*cfg.Req)-sinc0(q*cfg.Rgs));
c.reference=reference;
[cleanFull,cleanOsc]=signals(reference,c);
c.cleanOsc=movmean(cleanOsc,cfg.hodmd.pre.qSmoothWindow,1,'Endpoints','shrink');
c.cleanFull=movmean(cleanFull,cfg.hodmd.pre.qSmoothWindow,1,'Endpoints','shrink');
c.trueR=distance(reference,t,cfg.Req);
[~,SS,V]=svd(cleanFull,'econ');rv=SS*V';baseSigma=range(rv(1,:))/4;
generationCheck=table();
for n=1:numel(raw)
    delta=raw{n}.A-cleanFull;
    observedSigma=std(delta(:)); expectedSigma=baseSigma*cfg.noiseLevels(n)/100;
    meanDifference=norm(mean(raw{n}.A,3)-cleanFull,'fro')/norm(cleanFull,'fro');
    generationCheck=[generationCheck;table(cfg.noiseLevels(n),observedSigma,expectedSigma,meanDifference, ...
        'VariableNames',{'Noise_pct','ObservedSigma','ExpectedSigma','MeanSignalRelativeError'})]; %#ok<AGROW>
    assert(abs(observedSigma/expectedSigma-1)<0.05,'Generator/noise level mismatch.');
end
disp(generationCheck);
runDir=fullfile(cfg.outputRoot,[char(datetime('now','Format','yyyyMMdd_HHmmss_SSS')) '_' cfg.runLabel]);
assert(~isfolder(runDir),'Output already exists.');mkdir(runDir);
records=struct([]);rows=struct([]);
% Fit only the requested noisy input, without a noise-free control fit.
for n=1:numel(raw)
    seedList=cfg.seeds;noise=cfg.noiseLevels(n);
    for seed=seedList
        A=raw{n}.A(:,:,seed);
        timer=tic;
        h=hodmd_core('analyze',A,q,t,cfg.hodmd);
        targetSeconds=toc(timer);
        idx=h.selectedIndices;
        valid=~isempty(idx) && all(h.modeTable.DecayTime_fs(idx)>0) && all(isfinite(h.reconstruction(:)));
        target=h.reconstruction;
        fprintf('\nNoise %d%% seed %d, target mode ids %s\n',noise,seed,mat2str(idx));
        fits=cell(1,1);
        for method=1
            cc=c;cc.target=target;cc.residualScale=max(norm(target,'fro')/sqrt(numel(target)),eps);
            if valid
                positive=idx(h.modeTable.Frequency(idx)>0);
                assert(numel(positive)==1,'Exactly one oscillatory conjugate pair is required.');
                cc.fixedTime=[h.modeTable.Period_fs(positive),h.modeTable.DecayTime_fs(positive)];
                assert(all(isfinite(cc.fixedTime)) && all(cc.fixedTime>0),'Invalid fixed T/tau.');
                fprintf('Fixed HODMD T = %.6f fs, tau = %.6f fs\n',cc.fixedTime);
                fits{method}=fit_one(method,cc,starts);
            else
                fits{method}=struct('valid',false,'message','No decaying target pair in configured period window.');
            end
            row=summarize_fit(fits{method},method,noise,seed,c);
            if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<SAGROW>
        end
        rec=struct('noisePercent',noise,'seed',seed,'targetValid',valid,'experimentalHODMD',h, ...
            'targetSeconds',targetSeconds,'fits',{fits});
        if isempty(records),records=rec;else,records(end+1)=rec;end %#ok<SAGROW>
        save(fullfile(runDir,'checkpoint.mat'),'cfg','reference','q','t','records','rows','starts','generationCheck','audit','-v7.3');
    end
end
perFit=struct2table(rows);
fprintf('\nPer-fit results\n');disp(perFit);
results=struct('cfg',cfg,'q',q,'t_fs',t,'reference',reference,'parameterNames',{parameterNames}, ...
    'starts',starts,'generationCheck',generationCheck,'sourceAudit',audit, ...
    'records',records,'perFit',perFit, ...
    'trueR',c.trueR,'cleanOscillation',c.cleanOsc,'cleanFull',c.cleanFull, ...
    'notes','Direct fit: amplitude/phase free; T/tau fixed to target HODMD eigenvalue. Rgs, Req and scale fixed.');
if cfg.makeFigures,plot_direct(results,runDir);end
save(fullfile(runDir,'direct_fit_results.mat'),'-struct','results','-v7.3');
fprintf('\nSaved direct fit: %s\n',fullfile(runDir,'direct_fit_results.mat'));

function [full,osc]=signals(p,c)
    r=distance(p,c.tf,c.cfg.Req);
    fine=2*c.q.*c.f2.*(sinc0(c.q*r)-sinc0(c.q*c.cfg.Req));
    fullFine=fine+c.staticFine;
    oscFine=conv2(fine,c.kernel,'same');
    fullFine=conv2(fullFine,c.kernel,'same');
    osc=oscFine(:,c.sampleIdx);full=fullFine(:,c.sampleIdx);
end

function R=distance(p,t,Req)
    R=Req+p(1)*exp(-t/p(4)).*cos(2*pi*t/p(3)-p(2));
end

function y=sinc0(x)
    y=ones(size(x));idx=abs(x)>1e-12;y(idx)=sin(x(idx))./x(idx);
end

function [res,model]=objective(z,method,c)
    assert(method==1,'Only Direct fitting is enabled.');
    p=[c.cfg.lower+z.*(c.cfg.upper-c.cfg.lower),c.fixedTime];
    [~,osc]=signals(p,c);
    model=movmean(osc,c.cfg.hodmd.pre.qSmoothWindow,1,'Endpoints','shrink');
    res=(model(:)-c.target(:))/c.residualScale;
end

function fit=fit_one(method,c,starts)
    cfg=c.cfg;span=cfg.upper-cfg.lower;
    opts=optimoptions('lsqnonlin','Display','off','FiniteDifferenceType','central', ...
        'FiniteDifferenceStepSize',cfg.finiteDifferenceStep, ...
        'FunctionTolerance',cfg.functionTolerance,'StepTolerance',cfg.stepTolerance, ...
        'OptimalityTolerance',cfg.optimalityTolerance,'TypicalX',ones(1,2), ...
        'MaxIterations',cfg.maxIterations,'MaxFunctionEvaluations',cfg.maxFunctionEvaluations);
    attempts=repmat(struct('parameters',nan(1,4),'SSE',Inf,'exitflag',-999,'functionCount',0,'message',''),size(starts,1),1);
    best=Inf;timer=tic;fit=struct('valid',false);
    for k=1:size(starts,1)
        try
            fun=@(z)objective(z,method,c);
            [z,fval,~,flag,out]=lsqnonlin(fun,(starts(k,:)-cfg.lower)./span,zeros(1,2),ones(1,2),opts);
            p=[cfg.lower+z.*span,c.fixedTime];
            attempts(k)=struct('parameters',p,'SSE',fval*c.residualScale^2,'exitflag',flag, ...
                'functionCount',out.funcCount,'message',out.message);
            if isfinite(fval) && fval<best
                best=fval;fit.valid=true;fit.parameters=p;fit.exitflag=flag;fit.bestStart=k;
            end
        catch err
            attempts(k).message=err.message;
        end
        fprintf('  method %d start %d/%d SSE %.7g exit %d\n',method,k,size(starts,1),attempts(k).SSE,attempts(k).exitflag);
    end
    fit.seconds=toc(timer);fit.attempts=attempts;
    if ~fit.valid,fit.message='All starts failed.';return;end
    p=fit.parameters;
    fit.fixedTime_fs=c.fixedTime;
    fit.freeParameterIndices=cfg.freeParameterIndices;
    assert(isequal(p(3:4),c.fixedTime),'Fixed HODMD T/tau changed.');
    [~,fit.predictor]=objective((p(1:2)-cfg.lower)./span,method,c);
    [full,osc]=signals(p,c);
    fit.physicalOscillation=movmean(osc,cfg.hodmd.pre.qSmoothWindow,1,'Endpoints','shrink');
    fit.physicalFull=movmean(full,cfg.hodmd.pre.qSmoothWindow,1,'Endpoints','shrink');
    fit.R=distance(p,c.t,cfg.Req);
    fit.targetSSE=sum((fit.predictor(:)-c.target(:)).^2);
    fit.cleanOscRelativeError=norm(fit.physicalOscillation-c.cleanOsc,'fro')/norm(c.cleanOsc,'fro');
    fit.distanceRMSE=sqrt(mean((fit.R-c.trueR).^2));
    fit.phaseError=atan2(sin(p(2)-c.reference(2)),cos(p(2)-c.reference(2)));
    fit.nearBound=any((p(1:2)-cfg.lower)./span<0.001 | (cfg.upper-p(1:2))./span<0.001);
end

function row=summarize_fit(fit,method,noise,seed,c)
    labels={'Direct','TheoryHODMD'};
    row=struct('Noise_pct',noise,'Seed',seed,'Method',string(labels{method}), ...
        'Valid',fit.valid,'Amplitude_A',NaN,'Phase_rad',NaN,'Period_fs',NaN,'Decay_fs',NaN, ...
        'AmpError_A',NaN,'PhaseError_rad',NaN,'PeriodError_fs',NaN,'DecayError_fs',NaN, ...
        'DistanceRMSE_A',NaN,'CleanOscRelError',NaN,'TargetSSE',NaN,'Seconds',NaN,'Exitflag',-999,'NearBound',false);
    if ~fit.valid,return;end
    p=fit.parameters;e=p-c.reference;
    row.Amplitude_A=p(1);row.Phase_rad=p(2);row.Period_fs=p(3);row.Decay_fs=p(4);
    row.AmpError_A=e(1);row.PhaseError_rad=fit.phaseError;row.PeriodError_fs=e(3);row.DecayError_fs=e(4);
    row.DistanceRMSE_A=fit.distanceRMSE;row.CleanOscRelError=fit.cleanOscRelativeError;
    row.TargetSSE=fit.targetSSE;row.Seconds=fit.seconds;row.Exitflag=fit.exitflag;row.NearBound=fit.nearBound;
end


function plot_direct(r,runDir)
    rec=r.records(1);fit=rec.fits{1};
    assert(fit.valid,'Direct fit failed; inspect saved checkpoint.');
    target=rec.experimentalHODMD.reconstruction;
    f(1)=newfig('I2 90% N1 Direct fit');tl=tiledlayout(f(1),1,3);
    maps={target,fit.predictor,target-fit.predictor};
    titles={'Input HODMD component','Direct fit','Residual'};
    limit=max(abs([target(:);fit.predictor(:);maps{3}(:)]));
    for j=1:3
        ax=nexttile(tl);imagesc(ax,r.q,r.t_fs,maps{j}.');
        set(ax,'YDir','reverse');axis(ax,'square');clim(ax,[-limit limit]);
        colorbar(ax);title(ax,titles{j});xlabel(ax,'q (A^{-1})');ylabel(ax,'Time (fs)');
    end
    f(2)=newfig('I2 90% N1 distance');ax=axes(f(2));
    plot(ax,r.t_fs,r.trueR,'k-','LineWidth',1.5,'DisplayName','Actual');hold(ax,'on');
    plot(ax,r.t_fs,fit.R,'r-','LineWidth',1.3,'DisplayName','Fitted');
    yline(ax,r.cfg.Req,'k--','DisplayName','Req');xlim(ax,[300 2000]);
    xlabel(ax,'Time (fs)');ylabel(ax,'R(I-I) (A)');grid(ax,'on');legend(ax,'Location','best');
    drawnow;
    if r.cfg.exportFigures
        for j=1:numel(f),exportgraphics(f(j),fullfile(runDir,sprintf('figure_%02d.png',j)),'Resolution',130);end
    end
end

function f=newfig(name)
    f=figure();set(f,'WindowStyle','docked','Color','w','Name',name);
end


