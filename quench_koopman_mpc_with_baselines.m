function result = quench_koopman_mpc_with_baselines(options)
%quench_koopman_mpc_with_baselines Current manuscript Algorithms 1 and 2, self-contained.
% Run with no arguments for the full seven-controller comparison.
% Figures are enabled automatically in desktop MATLAB; batch mode saves data.
%   quench_koopman_mpc_with_baselines
%   result = quench_koopman_mpc_with_baselines('quick');   % short integration check
%   checks = quench_koopman_mpc_with_baselines('verify');  % equation/run-path checks only
% Options is also a struct overriding named cfg settings, e.g.:
%   result = quench_koopman_mpc_with_baselines(struct('makePlots',false,'saveFigures',false));
% Requires MATLAB and Optimization Toolbox. No helper files are required.
% Offline evidence uses UNCONSTRAINED virtual trajectories (Algorithm 1).
% Online uses cumulative EDMDc, frozen-window profile evidence, and then the
% CURRENT constrained control move (Algorithm 2; manuscript eqs. 23-24).
% Baselines use the fixed offline predictor; online uses only a declared
% coarse kinematic/input prior and newly measured transitions. Metrics use a
% common physical equilibrium input, separate from each internal reference.
% The supplied SelfMPC/COAt/DiffTune variants and their budgets are retained.
    if nargin == 0, options = struct(); end
%% Example 3: industrial-scale gasifier quench level, polynomial EDMDc-MPC.
% Run this SCRIPT directly (F5). Requires Optimization Toolbox only.
% Same seven tuning/controller adaptations and three figures as Example 2.
% Published industrial operating data: Yuan et al., Process Engineering,
% 24(3), 315-325, 2024, DOI 10.12034/j.issn.1009-606X.223196, Tables 1-2.
% Literature provides operating data and a CFD model, NOT these two ODEs.
% This is a reduced mass-balance/flow-loop simulation, not a field-validated model.
% Unpublished geometry, evaporation and flow-loop data are stated assumptions.
% x1=(h-hNom)/hScale; x2=hdot/rateScale; u=(mCmd-mNom)/mScale.
% h: metres above slag outlet; mCmd: nominal-head blackwater demand (kg/s).
% All MPC predictors are one fixed polynomial EDMDc model; no CasADi.
% u_ref and the steady level-rate state are solved ONLY from the identified model.
% Costs use level (m), level rate (m/s), and flow deviation u-u_ref (t/h).
% The true simulator and plots convert these coordinates to physical units.
% DiffTune uses true simulator Jacobians, as in the supplied VDP template.
% COAt and SelfMPC remain the adaptations in that template, not claims of
% exact reproductions of the respective original algorithms.
% Outputs are local to this script: results_quench_koopman and Figures.

% Reference information: at sample k, hold xref(:,k) and uref(:,k)
% over the prediction horizon. Future reference steps are not previewed.
% All controller rollouts, Bayesian windows and DiffTune use this policy.

clc;
close all
% Existing user figures are left open.

% This script uses MathWorks quadprog and optimoptions. MOSEK also ships a
% function named quadprog, but its options interface is different.
mathWorksQPFolder = fullfile(matlabroot, 'toolbox', 'optim', 'optim');
mathWorksQPFile = fullfile(mathWorksQPFolder, 'quadprog.m');
assert(isfile(mathWorksQPFile), 'Optimization Toolbox quadprog is required.');
if ~strcmpi(which('quadprog'), mathWorksQPFile)
    addpath(mathWorksQPFolder, '-begin');
end
assert(strcmpi(which('quadprog'), mathWorksQPFile), ...
    'A third-party quadprog is shadowing the Optimization Toolbox solver.');

%% 1. Experiment settings.
% [L] Literature data; [A] explicitly assumed; [D] derived.
cfg.pressure = 4.0e6;              % [L] Chamber pressure (Pa), Table 1.
cfg.wetSyngasFlow = 6851;          % [L] m^3/h at operating conditions, NOT Nm^3/h.
cfg.waterInletTPH = 139.7;         % [L] Water inlet, t/h, Table 1; effective sump-feed anchor.
cfg.rhoWater = 854.9;              % [L] kg/m^3, Table 2.
cfg.rhoSyngas = 20.6;             % [L] kg/m^3, Table 2; operating context only.
cfg.hNom = 3.6;                  % [L] Published operating level, metres.
cfg.effectiveArea = 10;           % [A] m^2; excludes internal equipment.
cfg.evaporationTPH = 0;          % [A] ADDITIONAL sump net loss; mean upstream cooling is excluded.
cfg.blackwaterTPH = cfg.waterInletTPH-cfg.evaporationTPH; % [D] Steady balance.
cfg.flowTau = 8;                % [A] Inner blackwater flow-loop time constant (s).
cfg.downstreamPressure = 3.8e6;   % [A] Effective high-pressure receiver, Pa.
% Model boundary is the liquid sump after upstream gas cooling. Mean gas
% saturation/evaporation is not computed a second time inside this reduction.
% This lumped hydraulic link stays pressurized; it is NOT a flashing-valve
% sizing model for direct discharge from 4 MPa into a low-pressure flash tank.
cfg.g = 9.81;
cfg.massNom = cfg.blackwaterTPH*1000/3600; % [D] kg/s.
cfg.massScale = 1/3.6;  % [D] kg/s per unit u; u is a t/h deviation.
cfg.rateScale = 1;           % [D] x2 is the physical level rate in m/s.
cfg.hScale = 1;                % [D] x1 is the physical level deviation in m.
cfg.hBaseReference = 3.9;         % [A] Level setpoint above the published low level.
cfg.hDomain = [3.0 5.0];          % [A] Simulation validity envelope, not SIS limits.
cfg.trainingLevelRange = [3.2 4.8];
cfg.dt = 4;                     % Sample time: 4 s; level process is slow.
cfg.rkSubsteps = 4;
cfg.Np = 60;                    % Prediction interval: 240 s; covers the inventory transient.
cfg.N = 1000;                    % Test duration: 1200 s (20 min).
cfg.uMax = 0.5*cfg.blackwaterTPH;                   % Demand is 69.85--209.55 t/h at nominal pressure head.
cfg.x0 = [(3.6-cfg.hNom)/cfg.hScale;0];
cfg.trainingRateLimit = 1.1*cfg.massScale*cfg.uMax ...
    /(cfg.rhoWater*cfg.effectiveArea*cfg.rateScale); % Cover commanded-flow range.
cfg.Q0 = eye(2); cfg.R0 = 1;        % Actual Q=q*Q0, R=r*R0.
cfg.qInitial = 1000; cfg.rInitial = 1;
cfg.useEDMD = true;                % Fixed polynomial EDMDc identification.
cfg.nTraining = 4000; cfg.seed = 7;
% The prediction horizon and data domain follow this slow inventory process.
cfg.polynomialDegree = 4;           % Total degree; try 2, 3, 4, or 5.
cfg.polynomialScale = [1.5;0.003];        % Scale higher-order monomials.
cfg.includeConstant = true;        % Include the degree-zero observable.


% Proposed offline / online Bayesian settings.
cfg.offlineEpisodes = 10;
cfg.offlineSteps = 10;              % 100 fixed offline states by default.
cfg.maxIter = 100;                  % Safety cap; convergence normally stops earlier.
cfg.tolerance = 1e-3;
cfg.settleSteps = 4;                % Extra real updates to show a short settled tail.
cfg.energyFloor = 1e-12;            % Hold weights when aggregate energy is too small.
cfg.window = 10;                     % Online window: 10 samples, 40 s.
cfg.enableOnlinePause = true;
cfg.onlinePauseError = 0.01;
cfg.onlineResumeError = 0.02;

% Grid-search baseline.
cfg.gridRatios = logspace(0, 8, 17);

% SelfMPC baseline. The original SelfMPC ML controller requires a stochastic
% prediction model. The present deterministic plant therefore uses a virtual
% isotropic state-disturbance covariance only to form the ML covariance gain.
% Its scalar variance mainly scales the objective; the propagated covariance
% structure is what changes the relative weighting across prediction steps.
cfg.selfMPCNoiseStd = 0.1;
cfg.selfMPCRidge = 1e-9;

% Shared external closed-loop objective used ONLY by COAt-MPC and DiffTune-MPC:
%   Jext = sum ||e||^2 / Ex_manual + sum (u-u_ref)^2 / Eu_manual.
% This gives equal dimensionless weight to state tracking and input energy.
cfg.externalStateWeight = 1;
cfg.externalInputWeight = 1;

% COAt-MPC / SAFEOPT-style settings. We tune theta=log10(q/r) on a 1-D grid.
cfg.coatIterations = 50;
cfg.coatGridSize = 121;
cfg.coatBeta = 2.0;                  % GP confidence multiplier.
cfg.coatSafetyMultiplier = 1.25;    % Jext must remain <= 1.25*initial Jext with high confidence.
cfg.coatLengthScale = 0.12;         % Matern-5/2 length scale on normalized theta in [0,1].
cfg.coatKernelScale = 0.50;          % GP signal std as a fraction of |initial score|.
cfg.coatNoiseScale = 0.01;           % GP observation std as a fraction of |initial score|.

% DiffTune-MPC settings. Because multiplying q and r by the same positive
% scalar does not change the MPC solution, we tune only theta=log(q/r), fix r=1,
% and set q=exp(theta). This removes the scale non-identifiability.
cfg.diffTuneIterations = 20;
cfg.diffTuneLearningRate = 0.20;
cfg.diffTuneTolerance = 1e-3;
cfg.diffTuneActiveTolerance = 1e-6;
cfg.diffTuneGradientClip = Inf;      % Set finite value if desired; Inf matches vanilla gradient descent.

% Reference and plotting.
cfg.referenceMode = 'steps';         % 'regulation' or 'steps'.
cfg.stepTimes = [400 800 1500 2000];                % Seconds; used only in 'steps' mode.
cfg.stepValues = [4.1 3.9 4.5 4.2]; % Physical level targets (m).
cfg.gainChangeTime = Inf;            % Optional actuator-gain change; default disabled.
cfg.gainAfterChange = 0.7;           % Multiplier; predictor remains fixed.
cfg.makePlots = true;               % false: numerical run only.
cfg.saveFigures = true;
cfg.plotUntil = 3000;                  % Trajectory plot limit (s); Inf shows the complete test.
cfg.lineWidth = 2.1;
cfg.referenceLineWidth = 1.4;
% Keep the previously requested line-width layering for coincident curves.
cfg.trajectoryLineWidths = [4.6 2.2 3.0 1.8 2.4 3.2 1.2];
cfg.trajectoryReferenceLineWidth = 0.9;

% Current-algorithm numerical settings (no offline samples in Algorithm 2).
cfg.caseName = 'quench';
cfg.onlineInitialInputGain = -8e-6; % Coarse sign/scale prior; not fitted.
cfg = current_options(cfg,options);
preflight = verify_current_algorithm(cfg);
if strcmpi(cfg.mode,'verify'), result = preflight; return; end
assert(exist('quadprog', 'file') == 2, 'Optimization Toolbox is required.');
assert(cfg.dt > 0 && cfg.uMax > 0 && cfg.R0 > 0);
assert(cfg.useEDMD, 'This experiment requires the fixed EDMDc predictor.');
assert(cfg.polynomialDegree>=1 && cfg.polynomialDegree==fix(cfg.polynomialDegree));
assert(isequal(size(cfg.polynomialScale),[2 1]) && all(cfg.polynomialScale>0));
assert(islogical(cfg.includeConstant) && isscalar(cfg.includeConstant));
assert(cfg.rkSubsteps >= 1 && cfg.rkSubsteps == fix(cfg.rkSubsteps));
assert(cfg.Np >= 1 && cfg.Np == fix(cfg.Np) && cfg.N >= 1 && cfg.N == fix(cfg.N));
assert(cfg.window >= 1 && cfg.window == fix(cfg.window));
assert(isscalar(cfg.enableOnlinePause) && islogical(cfg.enableOnlinePause));
assert(isfinite(cfg.onlinePauseError) && cfg.onlinePauseError >= 0 ...
    && isfinite(cfg.onlineResumeError) && cfg.onlineResumeError > cfg.onlinePauseError);
assert(cfg.maxIter >= 1 && cfg.maxIter == fix(cfg.maxIter));
assert(cfg.settleSteps >= 0 && cfg.settleSteps == fix(cfg.settleSteps) ...
    && cfg.settleSteps < cfg.maxIter);
assert(cfg.tolerance > 0 && cfg.energyFloor >= 0 && all(cfg.gridRatios > 0));
assert(isequal(size(cfg.x0), [2 1]) && isequal(size(cfg.Q0), [2 2]));
assert(norm(cfg.Q0-cfg.Q0', 'fro') < 1e-12 && all(eig(cfg.Q0) > 0));
assert(cfg.qInitial > 0 && cfg.rInitial > 0 && cfg.gainAfterChange > 0);
assert(cfg.offlineEpisodes >= 2 && cfg.offlineEpisodes == fix(cfg.offlineEpisodes) ...
    && cfg.offlineSteps >= 1 && cfg.offlineSteps == fix(cfg.offlineSteps));
assert(cfg.nTraining >= 3 && cfg.nTraining == fix(cfg.nTraining) && cfg.plotUntil > 0);
assert(cfg.lineWidth > 0 && cfg.referenceLineWidth > 0);
assert(numel(cfg.stepTimes) == numel(cfg.stepValues) && all(diff(cfg.stepTimes) > 0));
assert(any(strcmp(cfg.referenceMode, {'regulation', 'steps'})));
assert(cfg.coatIterations >= 1 && cfg.coatGridSize >= 11 && cfg.coatBeta > 0);
assert(cfg.coatSafetyMultiplier >= 1 && cfg.coatLengthScale > 0);
assert(cfg.diffTuneIterations >= 1 && cfg.diffTuneLearningRate > 0);
assert(cfg.selfMPCNoiseStd > 0 && cfg.selfMPCRidge >= 0);

assert(cfg.pressure>cfg.downstreamPressure && cfg.rhoWater>0);
assert(cfg.blackwaterTPH>0 && cfg.effectiveArea>0 && cfg.flowTau>0);
assert(cfg.massNom-cfg.massScale*cfg.uMax>0 && cfg.rateScale>0);
assert(cfg.hBaseReference>cfg.hDomain(1) && cfg.hBaseReference<cfg.hDomain(2));
assert(all(cfg.stepValues>cfg.hDomain(1) & cfg.stepValues<cfg.hDomain(2)));
fprintf('Industrial anchors: p=%.1f MPa, inlet=%.1f t/h, rho=%.1f kg/m^3, h=%.1f m.\n', ...
    cfg.pressure/1e6,cfg.waterInletTPH,cfg.rhoWater,cfg.hNom);
fprintf('ASSUMED: area=%.1f m^2, additional sump loss=%.1f t/h, flow-loop tau=%.1f s.\n', ...
    cfg.effectiveArea,cfg.evaporationTPH,cfg.flowTau);

folder = fileparts(mfilename('fullpath'));
outputFolder = fullfile(folder, 'results_quench_current');
if ~isempty(cfg.outputFolder), outputFolder = cfg.outputFolder; end
figureFolder = fullfile(outputFolder, 'Figures');
if ~isfolder(outputFolder), mkdir(outputFolder); end

%% 2. Identify a fixed polynomial EDMDc predictor z_next=A*z+B*u.
rng(cfg.seed,'twister');
trainingH = cfg.trainingLevelRange(1)+diff(cfg.trainingLevelRange)*rand(1,cfg.nTraining);
X = [(trainingH-cfg.hNom)/cfg.hScale;cfg.trainingRateLimit*(2*rand(1,cfg.nTraining)-1)];
trainingU = cfg.uMax*(2*rand(1,cfg.nTraining)-1);
Y = quench_step(X,trainingU,cfg,1);
dictionary = make_polynomial_dictionary(cfg.polynomialDegree, ...
    cfg.polynomialScale,cfg.includeConstant);
Z = polynomial_lift(X,dictionary);
Znext = polynomial_lift(Y,dictionary);
nZ = size(Z,1);
assert(cfg.nTraining>=nZ+1, 'Too few snapshots for the dictionary.');
% Input normalization improves conditioning without changing the fitted
% model: the reported B acts on the normalized command u.
ABscaled = Znext/[Z;trainingU/cfg.uMax];
A = ABscaled(:,1:nZ); B = ABscaled(:,nZ+1)/cfg.uMax;
if cfg.includeConstant
    % The constant observable is exactly invariant, not a dynamic state.
    A(3,:) = 0; A(3,3) = 1; B(3) = 0;
end
% The first two observables are the process coordinates in m and m/s.
C = [eye(2),zeros(2,nZ-2)];
trainingRMSE = sqrt(mean(sum((C*(A*Z+B*trainingU)-Y).^2,1)));
fprintf('Polynomial EDMDc: degree %d, %d observables, training norm RMSE %.6g\n', ...
    cfg.polynomialDegree,nZ,trainingRMSE);

%% 3. Projected physical-state prediction Xpred=O*Psi(x)+G*U.
% Cost uses process-state errors and u-u_ref, not every observable.
O = zeros(2*cfg.Np,nZ); G = zeros(2*cfg.Np,cfg.Np);
for j = 1:cfg.Np
    rows = 2*j-1:2*j;
    O(rows,:) = C*A^j;
    for ell = 1:j
        G(rows,ell) = C*A^(j-ell)*B;
    end
end
mdl = struct('A',A,'B',B,'C',C,'dictionary',dictionary, ...
    'O',O,'G',G,'Qbar',kron(eye(cfg.Np),cfg.Q0), ...
    'Rbar',cfg.R0*eye(cfg.Np),'uMax',cfg.uMax);
mdl.K = G'*mdl.Qbar*G;
mdl.GQ = G'*mdl.Qbar;
mdl.qpOptions = optimoptions('quadprog', 'Display', 'off', ...
    'Algorithm', 'interior-point-convex', 'OptimalityTolerance', 1e-9, ...
    'ConstraintTolerance', 1e-9, 'MaxIterations', 200);

% SelfMPC-style stacked prediction covariance from virtual process noise.
% First-order lift of virtual normalized-state noise at the origin.
noiseLift = polynomial_jacobian(zeros(2,1),dictionary);
D = zeros(2*cfg.Np, 2*cfg.Np);
for j = 1:cfg.Np
    rows = 2*j-1:2*j;
    for ell = 1:j
        cols = 2*ell-1:2*ell;
        D(rows, cols) = C*A^(j-ell)*noiseLift;
    end
end
SigmaSelf = cfg.selfMPCNoiseStd^2*(D*D') + cfg.selfMPCRidge*eye(2*cfg.Np);
SigmaSelf = (SigmaSelf+SigmaSelf')/2;
Lself = chol(SigmaSelf, 'lower');
mdl.selfW = Lself'\(Lself\eye(2*cfg.Np));

%% 4. Proposed offline Bayesian weight generator.
episodeStarts = [(linspace(3.5,4.3,cfg.offlineEpisodes)-cfg.hNom)/cfg.hScale; ...
    0.0006/cfg.rateScale*(-1).^(0:cfg.offlineEpisodes-1)];
offlineTarget = [(cfg.hBaseReference-cfg.hNom)/cfg.hScale;0];
[offlineInputReference,offlineReferenceInfo] = make_input_reference(mdl,offlineTarget,cfg);
offlineTarget = offlineReferenceInfo.stateTargets(:,1);
offlineRefStack = repmat(offlineTarget,cfg.Np,1);
M = cfg.offlineEpisodes*cfg.offlineSteps;
offlineX = zeros(2, M); sample = 0;

for episode = 1:cfg.offlineEpisodes
    x = episodeStarts(:, episode);
    for k = 1:cfg.offlineSteps
        sample = sample+1; offlineX(:, sample) = x;
        [U, ~] = mpc_batch(mdl, x, offlineRefStack, cfg.qInitial, cfg.rInitial, repmat(offlineInputReference,cfg.Np,1));
        x = surrogate_step(mdl,x,U(1));
    end
end
offlineXref = repmat(offlineRefStack,1,M);
offlineUref = repmat(offlineInputReference,cfg.Np,M);
[qOffline,rOffline,offlineTrace,offlineStatus,firstToleranceIteration,convergenceIteration,actualUpdates] = ...
    tune_offline_virtual(mdl,cfg,offlineX,offlineXref,offlineUref);

%% 5. Grid-search ratio on separate regulation starts (original baseline).
validationStarts = [( [3.5 3.8 4.3]-cfg.hNom)/cfg.hScale; [-0.0006 0 0.0006]/cfg.rateScale];
gridCost = zeros(numel(cfg.gridRatios), size(validationStarts, 2));
gridCfg = cfg; gridCfg.gainChangeTime = Inf;
zeroReference = repmat(offlineTarget,1,cfg.N+cfg.Np+1);
for j = 1:numel(cfg.gridRatios)
    for trial = 1:size(validationStarts, 2)
        run = simulate_mpc(mdl, gridCfg, [cfg.gridRatios(j) 1], false, ...
            validationStarts(:, trial), zeroReference);
        gridCost(j, trial) = run.Jx;
    end
end
[~, best] = min(mean(gridCost, 2)); gridRatio = cfg.gridRatios(best);
fprintf('Grid search selected q/r=%.9g\n', gridRatio);

%% 6. Common reference for all closed-loop comparisons.
xref = repmat([(cfg.hBaseReference-cfg.hNom)/cfg.hScale;0],1,cfg.N+cfg.Np+1);
if strcmp(cfg.referenceMode, 'steps')
    for j = 1:numel(cfg.stepTimes)
        first = max(1, floor(cfg.stepTimes(j)/cfg.dt)+1);
        if first <= size(xref, 2), xref(1, first:end) = (cfg.stepValues(j)-cfg.hNom)/cfg.hScale; end
    end
end

[~, referenceInfo, xref] = make_input_reference(mdl, xref, cfg);
uref = evaluation_input_reference(xref,cfg);
fprintf('Surrogate references [h(m), rate(m/s), flow_ref(t/h), physical residual]:\n');
disp([cfg.hNom+cfg.hScale*referenceInfo.stateTargets(1,:)', ...
    cfg.rateScale*referenceInfo.stateTargets(2,:)', ...
    (cfg.massNom+cfg.massScale*referenceInfo.inputTargets')*3.6,referenceInfo.physicalResidual']);

% Manual run is evaluated first because its squared energies normalize the
% external task objective used by COAt-MPC and DiffTune-MPC.
manualRun = simulate_mpc(mdl, cfg, [cfg.qInitial cfg.rInitial], false, cfg.x0, xref);
[~, manualEx2, manualEu2] = external_score(manualRun, xref, cfg, [1 1]);
scoreScale = [max(manualEx2, eps), max(manualEu2, eps)];
manualJext = external_score(manualRun, xref, cfg, scoreScale);
fprintf('External baseline score (manual) = %.6g\n', manualJext);

%% 7. COAt-MPC / SAFEOPT-style closed-loop tuning of q/r.
[coatRatio, coatTrace] = tune_coat_ratio(mdl, cfg, cfg.x0, xref, scoreScale);
fprintf('COAt-MPC selected q/r=%.9g after %d closed-loop evaluations.\n', ...
    coatRatio, size(coatTrace, 1));

%% 8. DiffTune-MPC closed-loop sensitivity tuning of q/r.
[diffRatio, diffTrace] = tune_difftune_ratio(mdl, cfg, cfg.x0, xref, scoreScale);
fprintf('DiffTune-MPC selected q/r=%.9g after %d gradient iterations.\n', ...
    diffRatio, size(diffTrace, 1));

%% 9. Seven controllers on the identical plant/reference.
names = {'Manual', 'Grid search', 'Offline Bayesian', 'Online Bayesian', ...
    'SelfMPC', 'COAt-MPC', 'DiffTune-MPC'};
runs = cell(1, numel(names));
runs{1} = manualRun;
runs{2} = simulate_mpc(mdl, cfg, [gridRatio 1], false, cfg.x0, xref);
runs{3} = simulate_mpc(mdl, cfg, [qOffline rOffline], false, cfg.x0, xref);
runs{4} = simulate_mpc(mdl, cfg, [cfg.qInitial cfg.rInitial], true, cfg.x0, xref);
runs{5} = simulate_selfmpc(mdl, cfg, cfg.x0, xref);
runs{6} = simulate_mpc(mdl, cfg, [coatRatio 1], false, cfg.x0, xref);
runs{7} = simulate_mpc(mdl, cfg, [diffRatio 1], false, cfg.x0, xref);

% Performance metrics and final/learned q/r where defined.
ratioValues = [cfg.qInitial/cfg.rInitial; gridRatio; qOffline/rOffline; ...
    runs{4}.q(end)/runs{4}.r(end); NaN; coatRatio; diffRatio];
values = zeros(numel(names), 10);
for j = 1:numel(names)
    s = runs{j};
    values(j, :) = [ratioValues(j) s.Jx s.Ju s.terminalError max(abs(s.u)) ...
        mean(abs(s.u) >= cfg.uMax-1e-7) s.meanStepSeconds s.heldUpdates ...
        nnz(s.paused) s.resumeCount];
end
metrics = array2table(values, 'VariableNames', {'WeightRatio', 'Jx', 'Ju', ...
    'TerminalError', 'PeakInput', 'SaturationFraction', 'MeanStepSeconds', ...
    'HeldUpdates', 'PausedUpdates', 'ResumeCount'});
metrics = addvars(metrics, names', 'Before', 1, 'NewVariableNames', 'Controller');
disp(metrics);

writetable(metrics, fullfile(outputFolder, 'metrics.csv'));
write_current_diagnostics(runs,names,cfg,outputFolder,mdl,offlineX,offlineXref,offlineUref,offlineTrace);
result = struct('cfg',cfg,'metrics',metrics,'runs',{runs},'names',{names}, ...
    'offlineTrace',offlineTrace,'offlineStatus',offlineStatus,'preflight',preflight, ...
    'outputFolder',outputFolder,'xref',xref,'evaluationUref',uref);

writetable(array2table([referenceInfo.stateTargets',referenceInfo.inputTargets', ...
    referenceInfo.equilibriumResidual',referenceInfo.liftMismatch'], ...
    'VariableNames',{'x1Ref','x2Ref','uRef','EquilibriumResidual','LiftMismatch'}), ...
    fullfile(outputFolder,'reference_targets.csv'));
writetable(array2table(offlineTrace,'VariableNames', ...
    {'Iteration','q','r','SurrogateLogEvidence','EX','EU','h','VirtualOutsideBoundsFraction','RelativeStep'}), ...
    fullfile(outputFolder,'offline_trace.csv'));
writematrix(coatTrace, fullfile(outputFolder, 'coat_trace.csv'));
writematrix(diffTrace, fullfile(outputFolder, 'difftune_trace.csv'));

% Human-readable records for the exact step experiment.
writetable(array2table([cfg.gridRatios(:),gridCost,mean(gridCost,2)], ...
    'VariableNames',{'WeightRatio','JxStart1','JxStart2','JxStart3','MeanJx'}), ...
    fullfile(outputFolder,'grid_search.csv'));
online = runs{4};
writetable(table((0:cfg.N-1)'*cfg.dt,online.q,online.r, ...
    online.q./online.r,online.paused,online.held, ...
    'VariableNames',{'Time','q','r','WeightRatio','Paused','Held'}), ...
    fullfile(outputFolder,'online_weights.csv'));
for j = 1:numel(names)
    s = runs{j};
    writematrix([(0:cfg.N)'*cfg.dt,s.x',[s.u NaN]',xref(:,1:cfg.N+1)',uref(:,1:cfg.N+1)'], ...
        fullfile(outputFolder,sprintf('trajectory_%d.csv',j)));
end
% Physical units for the process engineer and manuscript figures.
physicalMetrics=zeros(numel(names),6);
for j=1:numel(names)
    h=cfg.hNom+cfg.hScale*runs{j}.x(1,:);
    actualTPH=3.6*(cfg.massNom-cfg.rhoWater*cfg.effectiveArea*cfg.rateScale*runs{j}.x(2,:));
    commandTPH=3.6*(cfg.massNom+cfg.massScale*runs{j}.u);
    physicalMetrics(j,:)=[min(h),max(h),abs(h(end)-(cfg.hNom+cfg.hScale*xref(1,cfg.N+1))), ...
        min(actualTPH),max(actualTPH),max(commandTPH)];
    physical=table((0:cfg.N)'*cfg.dt,h',1000*cfg.rateScale*runs{j}.x(2,:)',actualTPH', ...
        [commandTPH NaN]',(cfg.hNom+cfg.hScale*xref(1,1:cfg.N+1))', ...
        (3.6*(cfg.massNom+cfg.massScale*uref(1:cfg.N+1)))', ...
        'VariableNames',{'TimeSeconds','LevelM','LevelRateMmPerS','ActualBlackwaterTPH', ...
        'DemandTPH','ReferenceLevelM','ReferenceDemandTPH'});
    writetable(physical,fullfile(outputFolder,sprintf('physical_trajectory_%d.csv',j)));
end
physicalMetricsTable=array2table(physicalMetrics,'VariableNames',{'MinimumLevelM','MaximumLevelM', ...
    'TerminalLevelErrorM','MinimumBlackwaterTPH','MaximumBlackwaterTPH','PeakDemandTPH'});
physicalMetricsTable=addvars(physicalMetricsTable,names','Before',1,'NewVariableNames','Controller');
writetable(physicalMetricsTable,fullfile(outputFolder,'physical_metrics.csv'));
writetable(table((cfg.hNom+cfg.hScale*referenceInfo.stateTargets(1,:))', ...
    (cfg.rateScale*referenceInfo.stateTargets(2,:))', ...
    (3.6*(cfg.massNom+cfg.massScale*referenceInfo.inputTargets))', ...
    referenceInfo.physicalResidual',referenceInfo.equilibriumResidual', ...
    referenceInfo.predictionResidual', ...
    'VariableNames',{'LevelM','LevelRateMPerS','ReferenceDemandTPH', ...
    'PhysicalResidual','LiftedResidual','HorizonResidual'}), ...
    fullfile(outputFolder,'physical_reference_targets.csv'));

writetable(table(cfg.dt,cfg.N,cfg.uMax,cfg.plotUntil,cfg.qInitial,cfg.rInitial, ...
    cfg.pressure,cfg.waterInletTPH,cfg.rhoWater,cfg.effectiveArea,cfg.evaporationTPH, ...
    cfg.flowTau,cfg.hNom,cfg.hScale,cfg.massNom,cfg.massScale,cfg.downstreamPressure, ...
    string(mat2str(cfg.stepTimes)),string(mat2str(cfg.stepValues)),cfg.polynomialDegree,cfg.Np, ...
    'VariableNames',{'dt','N','uMax','PlotUntil','qInitial','rInitial','PressurePa', ...
    'WaterInletTPH','WaterDensity','EffectiveArea','EvaporationTPH','FlowLoopTau', ...
    'hNom','hScale','massNom','massScale','DownstreamPressurePa','StepTimes','StepValues', ...
    'PolynomialDegree','Np'}),fullfile(outputFolder,'experiment_settings.csv'));

save(fullfile(outputFolder, 'quench_koopman_results.mat'), 'cfg', 'mdl', ...
    'offlineX', 'offlineTrace', 'offlineStatus', 'firstToleranceIteration', ...
    'convergenceIteration', 'actualUpdates', 'gridCost', 'gridRatio', ...
    'coatRatio', 'coatTrace', 'diffRatio', 'diffTrace', 'scoreScale', ...
    'runs', 'names', 'metrics', 'xref', 'uref', 'referenceInfo', 'X', 'Y', 'trainingU', ...
    'dictionary','C','trainingRMSE', '-v7');

%% 10. Three figures, matching Example 1.
if cfg.makePlots
%% 10. Figures.
set(groot, 'defaultAxesFontName', 'Times New Roman', 'defaultAxesFontSize', 10, ...
    'defaultAxesLineWidth', 0.8, 'defaultAxesBox', 'on', 'defaultFigureColor', 'w', ...
    'defaultTextInterpreter', 'latex', 'defaultLegendInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex');
legendStyle = {'Interpreter', 'latex', 'Box', 'off', 'FontSize', 10};
colors = [0 .4470 .7410; .8500 .3250 .0980; .9290 .6940 .1250; ...
    .4940 .1840 .5560; .4660 .6740 .1880; .3010 .7450 .9330; .6350 .0780 .1840];
styles = {'-', '--', '-.', ':', '-', '--', '-.'};
t = (0:cfg.N)*cfg.dt;
figures = gobjects(1, 3);

% (1) All closed-loop trajectories.
figures(1) = figure('Name', 'quench chamber polynomial EDMDc-MPC trajectories', 'Position', [100 100 1000 620]);
[~, trajectoryOrder] = sort(cfg.trajectoryLineWidths, 'descend');
for panel = 1:3
    subplot(3, 1, panel); hold on; grid on;
    trajectoryHandles = gobjects(1, numel(names));
    % Draw wider curves first; retain method order explicitly in the legend.
    for j = trajectoryOrder
        if panel <= 2
            trajectoryHandles(j) = plot(t, physical_state(runs{j}.x,cfg,panel), styles{j}, ...
                'Color', colors(j, :), 'LineWidth', cfg.trajectoryLineWidths(j));
        else
            trajectoryHandles(j) = plot(t, 3.6*(cfg.massNom+cfg.massScale*[runs{j}.u runs{j}.u(end)]), styles{j}, ...
                'Color', colors(j, :), 'LineWidth', cfg.trajectoryLineWidths(j));
        end
    end
    if panel <= 2
        referenceHandle = plot(t, physical_state(xref(:,1:cfg.N+1),cfg,panel), 'k--', ...
            'LineWidth', cfg.trajectoryReferenceLineWidth);
        if panel==1, ylabel('$h$ (m)'); else, ylabel('$\dot{m}_{\mathrm{out}}$ (t/h)'); end
    else
        plot(t, 3.6*(cfg.massNom+cfg.massScale*uref(:,1:cfg.N+1)), 'k--', 'LineWidth', cfg.trajectoryReferenceLineWidth);
        yline(3.6*(cfg.massNom+cfg.massScale*cfg.uMax), 'k:', 'LineWidth', cfg.trajectoryReferenceLineWidth);
        yline(3.6*(cfg.massNom-cfg.massScale*cfg.uMax), 'k:', 'LineWidth', cfg.trajectoryReferenceLineWidth);
        ylabel('$m_{\mathrm{cmd}}$ (t/h)', 'Interpreter', 'latex'); xlabel('Time (s)', 'Interpreter', 'latex');
    end
    xlim([t(1) min(t(end), cfg.plotUntil)]);
    if panel == 1
        legend([trajectoryHandles referenceHandle], [names {'Reference'}], ...
            'Location', 'northoutside', 'NumColumns', 4, ...
            legendStyle{:});
    end
end

% (2) Proposed offline evidence/fixed-point convergence.
figures(2) = figure('Name', 'Offline Bayesian weights', 'Position', [100 100 1000 300]);
subplot(1, 3, 1);
semilogy(offlineTrace(:, 1), offlineTrace(:, 2), '-o', ...
    'LineWidth', cfg.lineWidth, 'MarkerSize', 4); hold on;
semilogy(offlineTrace(:, 1), offlineTrace(:, 3), '--s', ...
    'LineWidth', cfg.lineWidth, 'MarkerSize', 4);
grid on; xlabel('Iteration', 'Interpreter', 'latex'); ylabel('Weight', 'Interpreter', 'latex'); legend({'$q$', '$r$'}, legendStyle{:});
xlim([0 max(1, actualUpdates)]);
subplot(1, 3, 2);
semilogy(offlineTrace(:, 1), offlineTrace(:, 2)./offlineTrace(:, 3), '-o', ...
    'LineWidth', cfg.lineWidth, 'MarkerSize', 4);
grid on; xlabel('Iteration', 'Interpreter', 'latex'); ylabel('$q/r$', 'Interpreter', 'latex');
if isfinite(convergenceIteration)
    xline(convergenceIteration, 'k:', 'Converged', 'LineWidth', cfg.referenceLineWidth, 'Interpreter', 'latex');
end
xlim([0 max(1, actualUpdates)]);
subplot(1, 3, 3);
plot(offlineTrace(:, 1), offlineTrace(:, 4)-offlineTrace(1, 4), '-o', ...
    'LineWidth', cfg.lineWidth, 'MarkerSize', 4);
grid on; xlabel('Iteration', 'Interpreter', 'latex'); ylabel('Shifted full-space surrogate log evidence', 'Interpreter', 'latex');
xlim([0 max(1, actualUpdates)]);

% (3) Online adaptation and closed-loop trade-off in one 1-by-3 figure.
figures(3) = figure('Name', 'Online adaptation and tracking versus effort', ...
    'Position', [100 100 1200 360]);
combinedLayout = tiledlayout(figures(3), 1, 3, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
weightTime = t(1:end-1);

% (a) The two online weights retain their original logarithmic scale.
weightAxes = nexttile(combinedLayout, 1);
semilogy(weightAxes, weightTime, runs{4}.q, '-', ...
    weightTime, runs{4}.r, '--', 'LineWidth', cfg.lineWidth);
grid(weightAxes, 'on');
xlim(weightAxes, [t(1) t(end)]);
xlabel(weightAxes, 'Time (s)', 'Interpreter', 'latex');
ylabel(weightAxes, 'Online weight', 'Interpreter', 'latex');
title(weightAxes, '(a) Online weights', 'Interpreter', 'latex');
legend(weightAxes, {'$q$', '$r$'}, 'Location', 'best', legendStyle{:});

% (b) Method colors match the points in panel (c).
ratioAxes = nexttile(combinedLayout, 2);
semilogy(ratioAxes, weightTime, runs{4}.q./runs{4}.r, '-', ...
    'Color', colors(4, :), 'LineWidth', cfg.lineWidth);
hold(ratioAxes, 'on');
yline(ratioAxes, qOffline/rOffline, '-.', 'Color', colors(3, :), ...
    'LineWidth', cfg.referenceLineWidth);
yline(ratioAxes, cfg.qInitial/cfg.rInitial, ':', 'Color', colors(1, :), ...
    'LineWidth', cfg.referenceLineWidth);
yline(ratioAxes, coatRatio, '--', 'Color', colors(6, :), ...
    'LineWidth', cfg.referenceLineWidth);
yline(ratioAxes, diffRatio, '-', 'Color', colors(7, :), ...
    'LineWidth', cfg.referenceLineWidth);
set(ratioAxes, 'YScale', 'log');
grid(ratioAxes, 'on');
xlim(ratioAxes, [t(1) t(end)]);
xlabel(ratioAxes, 'Time (s)', 'Interpreter', 'latex');
ylabel(ratioAxes, '$q/r$', 'Interpreter', 'latex');
title(ratioAxes, '(b) Weight ratio', 'Interpreter', 'latex');

% (c) Use the newly simulated common-reference Jx--Ju data for all seven methods.
tradeoffAxes = nexttile(combinedLayout, 3);
hold(tradeoffAxes, 'on');
grid(tradeoffAxes, 'on');
methodHandles = gobjects(1, numel(names));
markers = {'o','s','^','d','v','+','x'};
for j = 1:numel(names)
    methodHandles(j) = plot(tradeoffAxes, values(j, 2), values(j, 3), markers{j}, ...
        'Color', colors(j, :), 'MarkerFaceColor', colors(j, :), ...
        'MarkerSize', 8, 'LineWidth', cfg.lineWidth);
end
% Distinct method markers remain legible when configurations lie close together.
xlabel(tradeoffAxes, '$J_x$', 'Interpreter', 'latex');
ylabel(tradeoffAxes, '$J_u$', 'Interpreter', 'latex');
title(tradeoffAxes, '(c) Tracking and input deviation', 'Interpreter', 'latex');
methodLegend = legend(tradeoffAxes, methodHandles, names, ...
    'NumColumns', 4, legendStyle{:});
methodLegend.Layout.Tile = 'north';
% methodLegend.Title.String = 'Methods in panels (b) and (c)';
methodLegend.Title.Interpreter = 'latex';



if cfg.saveFigures
    if ~isfolder(figureFolder), mkdir(figureFolder); end
    stems = {'quench_trajectory','quench_offline_convergence', ...
        'quench_online_weights_tracking_vs_effort'};
    for j = 1:numel(figures)
        set(figures(j),'PaperPositionMode','auto');
        savefig(figures(j),fullfile(figureFolder,[stems{j} '.fig']));
        print(figures(j),fullfile(figureFolder,[stems{j} '.png']),'-dpng','-r300');
        print(figures(j),fullfile(figureFolder,[stems{j} '.eps']),'-depsc','-painters'); %#ok<PRTPT>
    end
    fprintf('Figures saved in %s\n',figureFolder);
end
end
fprintf('Results saved in %s\n',outputFolder);

%% Local function 1: batched MPC and Bayesian-update quantities.
end % Main comparison function.

function [U, st] = mpc_batch(mdl, X, Xref, q, r, Uref, applyBounds)
%mpc_batch Actual constrained input; Bayesian statistics ALWAYS use virtual input.
% Pass applyBounds=false for Algorithm 1, avoiding all constrained QP solves.
    if nargin < 6, Uref = zeros(size(mdl.G,2),size(X,2)); end
    if nargin < 7, applyBounds = true; end
    M = size(X,2); nU = size(mdl.G,2); nX = size(mdl.G,1);
    assert(isequal(size(Uref),[nU,M]) && all(isfinite(Uref(:))));
    assert(isfinite(q) && isfinite(r) && min(q,r) > 0);
    scale = max(q,r);
    H = (q/scale)*mdl.K+(r/scale)*mdl.Rbar; H = (H+H')/2;
    Z = current_lift(mdl,X);
    b = (q/scale)*mdl.GQ*(Xref-mdl.O*Z)+(r/scale)*mdl.Rbar*Uref;
    L = chol(H,'lower'); virtualU = L'\(L\b); U = virtualU;
    if applyBounds
        objectiveScale = max(norm(H,inf),realmin);
        for column = find(any(abs(U) > mdl.uMax+1e-9,1))
            [solution,~,flag] = quadprog(2*H/objectiveScale,-2*b(:,column)/objectiveScale, ...
                [],[],[],[],-mdl.uMax*ones(nU,1),mdl.uMax*ones(nU,1),[],mdl.qpOptions);
            assert(flag > 0 && all(isfinite(solution)),'Constrained MPC QP failed.');
            U(:,column) = solution;
        end
        boundViolation = max([0;abs(U(:))-mdl.uMax]);
        assert(all(isfinite(U(:))) && boundViolation <= 1e-7,'MPC input constraint violated.');
        U = min(max(U,-mdl.uMax),mdl.uMax);
    end
    if nargout > 1
        residual = Xref-mdl.O*Z-mdl.G*virtualU; inputResidual = virtualU-Uref;
        st.EX = sum(residual.*(mdl.Qbar*residual),'all');
        st.EU = sum(inputResidual.*(mdl.Rbar*inputResidual),'all');
        st.h = M/2*trace(L'\(L\mdl.Rbar))/scale;
        gamma2 = M*(q/scale)*trace(L'\(L\mdl.K));
        st.gamma = [M*(nX-nU)+2*r*st.h;gamma2];
        logdet2H = nU*log(2*scale)+2*sum(log(diag(L)));
        st.score = M*nX/2*log(q)+M*nU/2*log(r)-q*st.EX-r*st.EU-M/2*logdet2H;
        st.active = mean(abs(virtualU(:)) >= mdl.uMax-1e-7);
    end
end

% Current implementation is appended below.

function out = simulate_selfmpc(mdl, cfg, x0, xref)
    % SelfMPC has a predictive-density state cost and no separate input penalty.
    out.uref = evaluation_input_reference(xref(:,1:cfg.N),cfg);
    out.x = zeros(2, cfg.N+1); out.x(:,1) = x0;
    out.u = zeros(1, cfg.N); out.q = nan(cfg.N,1); out.r = nan(cfg.N,1);
    out.held = false(cfg.N,1); out.paused = false(cfg.N,1);
    out.trackingRMS = nan(cfg.N,1); out.resumeCount = 0; elapsed = 0; out.stepSeconds = zeros(cfg.N,1);
    H = mdl.G'*mdl.selfW*mdl.G; H = (H+H')/2;
    for k = 1:cfg.N
        timer = tic;
        % Only the currently effective target is available to the controller.
        future = repmat(xref(:, k), 1, cfg.Np);
        b = mdl.G'*mdl.selfW*(future(:)-mdl.O*polynomial_lift(out.x(:,k),mdl.dictionary));
        outU = solve_box_qp(H, b, mdl.uMax, mdl.qpOptions);
        out.u(k) = outU(1);
        stepElapsed = toc(timer); elapsed = elapsed+stepElapsed; out.stepSeconds(k) = stepElapsed;
        gain = 1;
        if (k-1)*cfg.dt >= cfg.gainChangeTime, gain = cfg.gainAfterChange; end
        out.x(:,k+1) = quench_step(out.x(:,k),out.u(k),cfg,gain);
        assert(all(isfinite(out.x(:,k+1))), 'Nonfinite SelfMPC plant state.');
    end
    errorX = out.x(:,1:cfg.N)-xref(:,1:cfg.N);
    out.Jx = sum(sqrt(sum(errorX.^2,1)));
    out.Ju = sum(abs(out.u-out.uref));
    out.terminalError = norm(out.x(:,end)-xref(:,cfg.N+1));
    out.heldUpdates = 0; out.meanStepSeconds = elapsed/cfg.N;
end

%% Local function 4: dimensionless external closed-loop task score.
function [score, Ex2, Eu2] = external_score(run, xref, cfg, scale)
    e = run.x(:,1:cfg.N)-xref(:,1:cfg.N);
    Ex2 = sum(e(:).^2);
    inputResidual = run.u-run.uref;
    Eu2 = sum(inputResidual(:).^2);
    score = cfg.externalStateWeight*Ex2/max(scale(1),eps) ...
        + cfg.externalInputWeight*Eu2/max(scale(2),eps);
end

%% Local function 5: COAt-MPC / SAFEOPT-style 1-D safe Bayesian optimization.
function [bestRatio, trace] = tune_coat_ratio(mdl, cfg, x0, xref, scoreScale)
    lo = log10(min(cfg.gridRatios)); hi = log10(max(cfg.gridRatios));
    thetaGrid = linspace(lo, hi, cfg.coatGridSize)';
    zGrid = (thetaGrid-lo)/max(hi-lo,eps);
    theta = min(max(log10(cfg.qInitial/cfg.rInitial), lo), hi);
    Xobs = zeros(cfg.coatIterations,1); Yobs = zeros(cfg.coatIterations,1);
    trace = nan(cfg.coatIterations,5); safetyCost = NaN; yPrior = NaN;
    for it = 1:cfg.coatIterations
        ratio = 10^theta;
        run = simulate_mpc(mdl, cfg, [ratio 1], false, x0, xref);
        score = external_score(run, xref, cfg, scoreScale);
        if it == 1
            safetyCost = cfg.coatSafetyMultiplier*score;
            yPrior = -score;
        end
        Xobs(it) = theta; Yobs(it) = -score;
        n = it;
        zObs = (Xobs(1:n)-lo)/max(hi-lo,eps);
        signalStd = cfg.coatKernelScale*max(1,abs(yPrior));
        noiseStd = cfg.coatNoiseScale*max(1,abs(yPrior));
        [mu, sd] = gp_matern52(zObs, Yobs(1:n), zGrid, yPrior, ...
            cfg.coatLengthScale, signalStd, noiseStd);
        lcb = mu-cfg.coatBeta*sd; ucb = mu+cfg.coatBeta*sd;
        safe = lcb >= -safetyCost;
        % Keep every actually observed safe point in the safe set.
        for jj = 1:n
            if -Yobs(jj) <= safetyCost
                [~, idxObs] = min(abs(thetaGrid-Xobs(jj)));
                safe(idxObs) = true;
            end
        end
        trace(it,:) = [it ratio score safetyCost interp1(thetaGrid,sd,theta,'linear','extrap')];
        if it == cfg.coatIterations, break; end
        if ~any(safe)
            [~, idxNext] = min(abs(thetaGrid-theta));
        else
            bestLCB = max(lcb(safe));
            potentialMax = safe & (ucb >= bestLCB);
            boundary = false(size(safe));
            ids = find(safe);
            for jj = ids(:)'
                if (jj > 1 && ~safe(jj-1)) || (jj < numel(safe) && ~safe(jj+1))
                    boundary(jj) = true;
                end
            end
            candidates = potentialMax | boundary;
            width = 2*cfg.coatBeta*sd;
            % Avoid exactly repeating an already evaluated grid point unless no alternative exists.
            sampled = false(size(thetaGrid));
            for jj = 1:n
                [~, kk] = min(abs(thetaGrid-Xobs(jj))); sampled(kk) = true;
            end
            available = candidates & ~sampled;
            if any(available)
                tmp = width; tmp(~available) = -Inf; [~,idxNext] = max(tmp);
            else
                available = safe & ~sampled;
                if any(available)
                    tmp = ucb; tmp(~available) = -Inf; [~,idxNext] = max(tmp);
                else
                    tmp = ucb; tmp(~safe) = -Inf; [~,idxNext] = max(tmp);
                end
            end
        end
        theta = thetaGrid(idxNext);
    end
    trace = trace(1:find(isfinite(trace(:,1)),1,'last'),:);
    safeObs = trace(:,3) <= trace(1,4);
    if any(safeObs)
        ids = find(safeObs); [~,jj] = min(trace(ids,3)); bestIdx = ids(jj);
    else
        [~,bestIdx] = min(trace(:,3));
    end
    bestRatio = trace(bestIdx,2);
end

%% Local function 6: Matern-5/2 GP posterior with constant prior mean.
function [mu, sd] = gp_matern52(x, y, xs, priorMean, ell, signalStd, noiseStd)
    x = x(:); y = y(:); xs = xs(:);
    D = abs(x-x');
    R = sqrt(5)*D/ell;
    K = signalStd^2*(1+R+R.^2/3).*exp(-R) + noiseStd^2*eye(numel(x));
    Ds = abs(x-xs'); Rs = sqrt(5)*Ds/ell;
    Ks = signalStd^2*(1+Rs+Rs.^2/3).*exp(-Rs);
    L = chol((K+K')/2 + 1e-12*eye(size(K)), 'lower');
    alpha = L'\(L\(y-priorMean));
    mu = priorMean + Ks'*alpha;
    v = L\Ks;
    varPost = max(signalStd^2-sum(v.^2,1)', 1e-12);
    sd = sqrt(varPost);
end

%% Local function 7: DiffTune-MPC closed-loop sensitivity descent.
function [ratio, trace] = tune_difftune_ratio(mdl, cfg, x0, xref, scoreScale)
    uref = make_input_reference(mdl,xref,cfg);
    evaluationUref = evaluation_input_reference(xref,cfg);
    thetaMin = log(min(cfg.gridRatios)); thetaMax = log(max(cfg.gridRatios));
    theta = min(max(log(cfg.qInitial/cfg.rInitial),thetaMin),thetaMax);
    trace = nan(cfg.diffTuneIterations,4);
    for it = 1:cfg.diffTuneIterations
        q = exp(theta); r = 1;
        x = x0; dx = zeros(2,1); loss = 0; grad = 0;
        for k = 1:cfg.N
            % Match the reference information used in closed-loop evaluation.
            future = repmat(xref(:,k), 1, cfg.Np);
            [U, duTheta, duDx] = mpc_first_sensitivity(mdl, x, future(:), q, r, ...
                cfg.diffTuneActiveTolerance, repmat(uref(:,k), cfg.Np, 1));
            u = U(1);
            duTotal = duTheta + duDx*dx;
            e = x-xref(:,k); eu = u-evaluationUref(k);
            loss = loss + cfg.externalStateWeight*(e'*e)/scoreScale(1) ...
                + cfg.externalInputWeight*(eu^2)/scoreScale(2);
            grad = grad + cfg.externalStateWeight*(2*e'*dx)/scoreScale(1) ...
                + cfg.externalInputWeight*(2*eu*duTotal)/scoreScale(2);
            gain = 1;
            if (k-1)*cfg.dt >= cfg.gainChangeTime, gain = cfg.gainAfterChange; end
            [xNext,Fx,Fu] = quench_step(x,u,cfg,gain);
            dxNext = Fx*dx+Fu*duTotal;
            assert(all(isfinite([xNext;dxNext])), 'Nonfinite quench chamber sensitivity.');
            x = xNext; dx = dxNext;
        end
        if isfinite(cfg.diffTuneGradientClip)
            grad = max(min(grad,cfg.diffTuneGradientClip),-cfg.diffTuneGradientClip);
        end
        trace(it,:) = [it exp(theta) loss grad];
        thetaNext = min(max(theta-cfg.diffTuneLearningRate*grad,thetaMin),thetaMax);
        if abs(thetaNext-theta) <= cfg.diffTuneTolerance
            theta = thetaNext;
            trace = trace(1:it,:);
            break
        end
        theta = thetaNext;
    end
    if any(~isfinite(trace(:,1)))
        trace = trace(1:find(isfinite(trace(:,1)),1,'last'),:);
    end
    ratio = exp(theta);
end

%% Local function 8: first-action sensitivities from the box-QP KKT system.
function [U, duTheta, duDx] = mpc_first_sensitivity(mdl, x, xref, q, r, activeTol, Uref)
    if nargin < 7, Uref = zeros(size(mdl.G,2),1); end
    H = q*mdl.K+r*mdl.Rbar; H = (H+H')/2;
    c = mdl.GQ*(xref-mdl.O*polynomial_lift(x,mdl.dictionary)); b = q*c+r*mdl.Rbar*Uref;
    U = solve_box_qp(H, b, mdl.uMax, mdl.qpOptions);
    active = abs(abs(U)-mdl.uMax) <= activeTol;
    free = ~active;
    dUtheta = zeros(size(U)); dUdx = zeros(numel(U),numel(x));
    if any(free)
        F = find(free); Aidx = find(active);
        Hff = H(F,F);
        dH = q*mdl.K; dbTheta = q*c;
        rhsTheta = dbTheta(F)-dH(F,F)*U(F);
        if ~isempty(Aidx), rhsTheta = rhsTheta-dH(F,Aidx)*U(Aidx); end
        dUtheta(F) = Hff\rhsTheta;
        dbdx = -q*mdl.GQ*mdl.O*polynomial_jacobian(x,mdl.dictionary);
        dUdx(F,:) = Hff\dbdx(F,:);
    end
    duTheta = dUtheta(1); duDx = dUdx(1,:);
end

%% Local function 9: positive-definite box-constrained QP solver.
function U = solve_box_qp(H, b, uMax, options)
    H = (H+H')/2;
    n = size(H,1);
    L = chol(H + 1e-14*eye(n), 'lower');
    U = L'\(L\b);
    if any(abs(U) > uMax+1e-9)
        objectiveScale = max(norm(H,inf),realmin);
        [U,~,flag] = quadprog(2*H/objectiveScale,-2*b/objectiveScale,[],[],[],[], ...
            -uMax*ones(n,1),uMax*ones(n,1),[],options);
        assert(flag > 0 && all(isfinite(U)), 'Box QP failed.');
    end
    boundViolation = max([0; abs(U(:))-uMax]);
    assert(all(isfinite(U(:))) && boundViolation <= 1e-7, ...
        'Box QP input bound violated by %.3g (allowed tolerance: 1e-7).', ...
        boundViolation);
    U = min(max(U, -uMax), uMax);
end


function [Y,Fx,Fu] = quench_step(X,U,cfg,gain)
    % Zero-order-hold input and RK4, vectorized over independent samples.
    h = cfg.dt/cfg.rkSubsteps; Y = X;
    needJacobian = nargout>1;
    if needJacobian
        assert(size(X,2)==1 && isscalar(U), 'Jacobian requires a single state.');
        Fx = eye(2); Fu = zeros(2,1);
    end
    for substep = 1:cfg.rkSubsteps
        k1 = quench_rhs(Y,U,cfg,gain);
        y2 = Y+(h/2)*k1; k2 = quench_rhs(y2,U,cfg,gain);
        y3 = Y+(h/2)*k2; k3 = quench_rhs(y3,U,cfg,gain);
        y4 = Y+h*k3; k4 = quench_rhs(y4,U,cfg,gain);
        if needJacobian
            level=cfg.hNom+cfg.hScale*Y(1);
            pressureRatio=(cfg.pressure-cfg.downstreamPressure+cfg.rhoWater*cfg.g*level) ...
                /(cfg.pressure-cfg.downstreamPressure+cfg.rhoWater*cfg.g*cfg.hNom);
            % Input derivative is state-dependent; evaluate at each RK4 stage.
            b = [0;-cfg.massScale*gain*sqrt(pressureRatio) ...
                /(cfg.rhoWater*cfg.effectiveArea*cfg.flowTau*cfg.rateScale)];
            J1 = quench_rhs_jacobian(Y,U,cfg,gain);
            J2 = quench_rhs_jacobian(y2,U,cfg,gain);
            J3 = quench_rhs_jacobian(y3,U,cfg,gain);
            J4 = quench_rhs_jacobian(y4,U,cfg,gain);
            k1x = J1; k1u = b;
            k2x = J2*(eye(2)+(h/2)*k1x); k2u = J2*((h/2)*k1u)+quench_input_jacobian(y2,cfg,gain);
            k3x = J3*(eye(2)+(h/2)*k2x); k3u = J3*((h/2)*k2u)+quench_input_jacobian(y3,cfg,gain);
            k4x = J4*(eye(2)+h*k3x); k4u = J4*(h*k3u)+quench_input_jacobian(y4,cfg,gain);
            subFx = eye(2)+(h/6)*(k1x+2*k2x+2*k3x+k4x);
            subFu = (h/6)*(k1u+2*k2u+2*k3u+k4u);
            Fu = subFx*Fu+subFu; Fx = subFx*Fx;
        end
        Y = Y+(h/6)*(k1+2*k2+2*k3+k4);
    end
    level=cfg.hNom+cfg.hScale*Y(1,:);
    massOut=(cfg.waterInletTPH-cfg.evaporationTPH)*1000/3600 ...
        -cfg.rhoWater*cfg.effectiveArea*cfg.rateScale*Y(2,:);
    assert(all(level>cfg.hDomain(1) & level<cfg.hDomain(2)), ...
        'Level left the stated model-validity domain; revise the test or model.');
    assert(all(massOut>0),'Unphysical negative blackwater discharge.');
end

function F = quench_rhs(X,U,cfg,gain)
    % Physical ODEs before coordinate conversion:
    % rho*A*hdot=m_in-m_evap-m_out;
    % tau*m_out_dot=(m_nom+m_scale*u)*sqrt(head(h)/head(hNom))-m_out.
    % m_nom=m_in-m_evap. Inner flow dynamics and head law are assumptions.
    h=cfg.hNom+cfg.hScale*X(1,:);
    nominalHead=cfg.pressure-cfg.downstreamPressure+cfg.rhoWater*cfg.g*cfg.hNom;
    pressureRatio=(cfg.pressure-cfg.downstreamPressure+cfg.rhoWater*cfg.g*h)/nominalHead;
    assert(all(pressureRatio>0),'Invalid hydraulic pressure difference.');
    demand=cfg.massNom+cfg.massScale*gain*U;
    netFeed=(cfg.waterInletTPH-cfg.evaporationTPH)*1000/3600;
    F=[cfg.rateScale/cfg.hScale*X(2,:); ...
        (netFeed-demand.*sqrt(pressureRatio) ...
         -cfg.rhoWater*cfg.effectiveArea*cfg.rateScale*X(2,:)) ...
        /(cfg.rhoWater*cfg.effectiveArea*cfg.flowTau*cfg.rateScale)];
end

function J = quench_rhs_jacobian(x,u,cfg,gain)
    h=cfg.hNom+cfg.hScale*x(1);
    nominalHead=cfg.pressure-cfg.downstreamPressure+cfg.rhoWater*cfg.g*cfg.hNom;
    head=cfg.pressure-cfg.downstreamPressure+cfg.rhoWater*cfg.g*h;
    demand=cfg.massNom+cfg.massScale*gain*u;
    J=[0,cfg.rateScale/cfg.hScale; ...
        -demand*cfg.g*cfg.hScale/(2*cfg.effectiveArea*cfg.flowTau*cfg.rateScale ...
        *nominalHead*sqrt(head/nominalHead)),-1/cfg.flowTau];
end

function b = quench_input_jacobian(x,cfg,gain)
    h=cfg.hNom+cfg.hScale*x(1);
    nominalHead=cfg.pressure-cfg.downstreamPressure+cfg.rhoWater*cfg.g*cfg.hNom;
    head=cfg.pressure-cfg.downstreamPressure+cfg.rhoWater*cfg.g*h;
    b=[0;-cfg.massScale*gain*sqrt(head/nominalHead) ...
        /(cfg.rhoWater*cfg.effectiveArea*cfg.flowTau*cfg.rateScale)];
end

function y = physical_state(x,cfg,panel)
    if panel == 1
        y = cfg.hNom+cfg.hScale*x(1,:);
    else
        y = 3.6*((cfg.waterInletTPH-cfg.evaporationTPH)*1000/3600 ...
            -cfg.rhoWater*cfg.effectiveArea*cfg.rateScale*x(2,:));
    end
end

%% Polynomial dictionary: raw x1,x2, optional 1, scaled degree-2..d terms.
function dictionary = make_polynomial_dictionary(degree,stateScale,includeConstant)
    nTerms = (degree+1)*(degree+2)/2-double(~includeConstant);
    powers = zeros(nTerms,2);
    powers(1:2,:) = [1 0;0 1]; count = 2;
    if includeConstant
        count=count+1; powers(count,:)=[0 0];
    end
    for totalDegree=2:degree
        for firstPower=totalDegree:-1:0
            count=count+1;
            powers(count,:)=[firstPower,totalDegree-firstPower];
        end
    end
    assert(count==nTerms);
    denominator = stateScale(1).^powers(:,1).*stateScale(2).^powers(:,2);
    denominator(1:2)=1;              % Physical state coordinates stay raw.
    dictionary=struct('degree',degree,'stateScale',stateScale, ...
        'includeConstant',includeConstant,'powers',powers,'denominator',denominator);
end

function Z = polynomial_lift(X,dictionary)
    assert(size(X,1)==2,'Expected physical states of size 2-by-M.');
    Z=(X(1,:).^dictionary.powers(:,1)).*(X(2,:).^dictionary.powers(:,2)) ...
        ./dictionary.denominator;
    assert(all(isfinite(Z(:))),'Nonfinite polynomial observables.');
end

function J = polynomial_jacobian(x,dictionary)
    assert(isequal(size(x),[2 1]),'A single physical state is required.');
    powers=dictionary.powers; J=zeros(size(powers,1),2);
    first=powers(:,1)>0; second=powers(:,2)>0;
    J(first,1)=powers(first,1).*x(1).^(powers(first,1)-1) ...
        .*x(2).^powers(first,2)./dictionary.denominator(first);
    J(second,2)=powers(second,2).*x(1).^powers(second,1) ...
        .*x(2).^(powers(second,2)-1)./dictionary.denominator(second);
end


%% Consistent steady target from the identified surrogate's physical rows.
% Current reference implementation is appended below.

function cfg = current_options(cfg, options)
%current_options Explicit implementation settings for Algorithms 1 and 2.
    cfg.mode = 'full';
    cfg.makePlots = usejava('desktop'); % Interactive MATLAB: plots; batch: numerical results.
    cfg.outputFolder = '';
    cfg.ratioBounds = [1, 1e10];
    cfg.ratioGridSize = 65;
    cfg.ratioLogTolerance = 1e-6;
    cfg.onlinePriorDelta = 1e6;
    cfg.onlineA0 = [];
    cfg.onlineB0 = [];
    cfg.onlineSensitivityFloor = 1e-14;
    cfg.onlinePredictionLimit = 1e8;
    cfg.onlineConditionFloor = 1e-14;
    if ischar(options) || isstring(options), options = struct('mode', char(options)); end
    assert(isstruct(options) && isscalar(options), 'Options must be a scalar struct or a mode string.');
    if isfield(options, 'mode') && strcmpi(options.mode, 'quick')
        cfg.N = 100;
        cfg.gridRatios = logspace(0, 8, 5);
        cfg.coatIterations = 3;
        cfg.diffTuneIterations = 3;
        cfg.offlineEpisodes = 4;
        cfg.offlineSteps = 5;
        cfg.makePlots = false;
        cfg.saveFigures = false;
        cfg.stepTimes = cfg.dt*floor(cfg.N*[0.3, 0.65]);
        cfg.stepValues = [cfg.stepValues(1), cfg.stepValues(end)];
    end
    fields = fieldnames(options);
    for k = 1:numel(fields)
        assert(isfield(cfg, fields{k}), 'Unknown option: %s', fields{k});
        cfg.(fields{k}) = options.(fields{k});
    end
    assert(numel(cfg.ratioBounds) == 2 && all(isfinite(cfg.ratioBounds)) ...
        && 0 < cfg.ratioBounds(1) && cfg.ratioBounds(1) < cfg.ratioBounds(2));
    assert(cfg.ratioGridSize >= 5 && cfg.ratioGridSize == fix(cfg.ratioGridSize));
    assert(isfinite(cfg.onlinePriorDelta) && cfg.onlinePriorDelta > 0);
    assert(any(strcmpi(cfg.mode, {'full', 'quick', 'verify'})), 'Unknown run mode.');
end

function Z = current_lift(mdl, X)
%current_lift Use only the fixed dictionary, without learned-model information.
    if isfield(mdl, 'dictionary'), Z = polynomial_lift(X, mdl.dictionary); else, Z = X; end
end

function xNext = surrogate_step(mdl, x, u)
%surrogate_step Generate offline states using the frozen surrogate, Algorithm 1.
    if isfield(mdl, 'C'), C = mdl.C; else, C = eye(2); end
    xNext = C*(mdl.A*current_lift(mdl, x)+mdl.B*u);
    assert(all(isfinite(xNext)), 'Nonfinite offline surrogate state.');
end

function [q, r, history, status, firstTolerance, convergence, updates] = ...
        tune_offline_virtual(mdl, cfg, X, Xref, Uref)
%tune_offline_virtual Algorithm 1: only unconstrained virtual residuals are used.
    q = cfg.qInitial; r = cfg.rInitial;
    history = nan(cfg.maxIter+1, 9);
    status = 'iteration limit'; firstTolerance = NaN; convergence = NaN;
    consecutive = 0; updates = 0;
    for iteration = 1:cfg.maxIter
        [~, stats] = mpc_batch(mdl, X, Xref, q, r, Uref, false);
        history(iteration,:) = [updates, q, r, stats.score, stats.EX, stats.EU, stats.h, stats.active, NaN];
        if min(stats.EX, stats.EU) <= cfg.energyFloor
            status = 'near-zero virtual energy: weights held'; break
        end
        proposal = stats.gamma./(2*[stats.EX; stats.EU]);
        if any(~isfinite(proposal)) || any(proposal <= 0)
            status = 'invalid virtual update: weights held'; break
        end
        relative = max(abs(proposal-[q;r])./[q;r]);
        history(iteration,9) = relative;
        q = proposal(1); r = proposal(2); updates = updates+1;
        if relative <= cfg.tolerance
            if isnan(firstTolerance), firstTolerance = updates; end
            consecutive = consecutive+1;
        else
            consecutive = 0;
        end
        if consecutive >= cfg.settleSteps+1
            convergence = updates-cfg.settleSteps; status = 'converged'; break
        end
    end
    [~, stats] = mpc_batch(mdl, X, Xref, q, r, Uref, false);
    history(updates+1,:) = [updates, q, r, stats.score, stats.EX, stats.EU, stats.h, stats.active, NaN];
    history = history(1:updates+1,:);
    fprintf('Offline Algorithm 1: %s; %d updates; q/r=%.9g\n', status, updates, q/r);
end

function model = prediction_matrices(model, cfg)
%prediction_matrices Project lifted predictions to the physical state space.
    nZ = size(model.A,1); nU = cfg.Np;
    model.O = zeros(2*cfg.Np,nZ); model.G = zeros(2*cfg.Np,nU);
    powers = cell(1,cfg.Np+1); powers{1} = eye(nZ);
    for j = 1:cfg.Np, powers{j+1} = powers{j}*model.A; end
    for j = 1:cfg.Np
        rows = 2*j-1:2*j;
        model.O(rows,:) = model.C*powers{j+1};
        for ell = 1:j, model.G(rows,ell) = model.C*powers{j-ell+1}*model.B; end
    end
    model.K = model.G'*model.Qbar*model.G; model.K = (model.K+model.K')/2;
    model.GQ = model.G'*model.Qbar;
end

function [model, identification] = initialize_online_model(base, cfg)
%initialize_online_model Cold initialization without offline transition samples.
% The declared nonzero input seed encodes a coarse sign/scale prior, not EDMDc data.
    nZ = size(base.A,1);
    A0 = eye(nZ); A0(1,2) = cfg.dt;
    B0 = zeros(nZ,1); B0(2) = cfg.onlineInitialInputGain;
    if ~isempty(cfg.onlineA0), A0 = cfg.onlineA0; end
    if ~isempty(cfg.onlineB0), B0 = cfg.onlineB0; end
    assert(isequal(size(A0),[nZ,nZ]) && isequal(size(B0),[nZ,1]));
    assert(all(isfinite(A0(:))) && all(isfinite(B0(:))));
    model = struct('A',A0,'B',B0,'C',[eye(2),zeros(2,nZ-2)], ...
        'Qbar',base.Qbar,'Rbar',base.Rbar,'uMax',base.uMax,'qpOptions',base.qpOptions);
    constantIndex = [];
    if isfield(base,'dictionary')
        model.dictionary = base.dictionary;
        if base.dictionary.includeConstant, constantIndex = 3; end
    end
    if ~isempty(constantIndex)
        model.A(constantIndex,:) = 0; model.A(constantIndex,constantIndex) = 1; model.B(constantIndex) = 0;
    end
    model = prediction_matrices(model,cfg);
    W0 = [model.A,model.B];
    identification = struct('P',cfg.onlinePriorDelta*eye(nZ+1), ...
        'information',eye(nZ+1)/cfg.onlinePriorDelta,'accumulator',W0/cfg.onlinePriorDelta, ...
        'W0',W0,'constantIndex',constantIndex,'samples',0);
end

function [candidate, identification, accepted, fitError] = ...
        recursive_online_update(model, identification, previousX, previousU, currentX, cfg)
%recursive_online_update Cumulative EDMDc update; all transitions are retained.
    z = current_lift(model,previousX); zNext = current_lift(model,currentX);
    phi = [z;previousU];
    fitError = norm(model.C*(model.A*z+model.B*previousU)-currentX);
    v = identification.P*phi; beta = 1+phi'*v;
    identification.P = identification.P-(v*v')/beta;
    identification.P = (identification.P+identification.P')/2;
    identification.information = identification.information+phi*phi';
    identification.accumulator = identification.accumulator+zNext*phi';
    identification.samples = identification.samples+1;
    % Solving the accumulated normal equations evaluates G_k P_k without
    % multiplying an ill-conditioned inverse. P_k is recomputed consistently.
    [L, flag] = chol((identification.information+identification.information')/2,'lower');
    candidate = model; accepted = false;
    if flag ~= 0 || ~isfinite(beta) || beta <= 0, return; end
    identification.P = L'\(L\eye(size(L)));
    W = (L'\(L\identification.accumulator'))';
    if ~isempty(identification.constantIndex)
        index = identification.constantIndex;
        W(index,:) = 0; W(index,index) = 1;
    end
    if any(~isfinite(W(:))), return; end
    proposed = model; proposed.A = W(:,1:end-1); proposed.B = W(:,end);
    proposed = prediction_matrices(proposed,cfg);
    if any(~isfinite(proposed.O(:))) || any(~isfinite(proposed.G(:))), return; end
    predictionMagnitude = max([abs(proposed.O(:));cfg.uMax*abs(proposed.G(:))]);
    if predictionMagnitude > cfg.onlinePredictionLimit, return; end
    H = proposed.K+proposed.Rbar/cfg.ratioBounds(2);
    if rcond(H) < cfg.onlineConditionFloor, return; end
    candidate = proposed; accepted = true;
end

function [uref, info, unchangedXref] = make_input_reference(mdl, xref, cfg)
%make_input_reference Algorithm 2: bounded physical one-step SISO target.
% No high-order observable or physical state reference is optimized here.
    if isfield(mdl,'C'), C = mdl.C; else, C = eye(2); end
    [targets,~,targetIndex] = unique(xref','rows','stable');
    n = size(targets,1); inputs = zeros(1,n);
    physicalResidual = zeros(1,n); liftedResidual = zeros(1,n); horizonResidual = zeros(1,n);
    degenerate = false(1,n); sensitivity = C*mdl.B;
    denominator = sensitivity'*cfg.Q0*sensitivity;
    for j = 1:n
        xr = targets(j,:)'; zr = current_lift(mdl,xr);
        offset = xr-C*mdl.A*zr;
        if isfinite(denominator) && denominator > cfg.onlineSensitivityFloor
            inputs(j) = min(max((sensitivity'*cfg.Q0*offset)/denominator,-cfg.uMax),cfg.uMax);
        else
            inputs(j) = 0; degenerate(j) = true; % Declared feasible minimum-norm fallback.
        end
        physicalResidual(j) = norm(C*(mdl.A*zr+mdl.B*inputs(j))-xr);
        liftedResidual(j) = norm((eye(size(mdl.A))-mdl.A)*zr-mdl.B*inputs(j));
        horizonResidual(j) = norm(mdl.O*zr+mdl.G*repmat(inputs(j),cfg.Np,1)-repmat(xr,cfg.Np,1));
    end
    uref = reshape(inputs(targetIndex),1,[]); unchangedXref = xref;
    info = struct('stateTargets',targets','inputTargets',inputs,'physicalResidual',physicalResidual, ...
        'equilibriumResidual',liftedResidual,'predictionResidual',horizonResidual, ...
        'liftMismatch',zeros(1,n),'degenerate',degenerate,'source','bounded physical one-step least squares');
end

function profile = prepare_profile(mdl, X, Xref, Uref)
%prepare_profile Frozen model and STORED references; stable output-space energy.
    SQ = chol(mdl.Qbar); SR = chol(mdl.Rbar);
    B = SQ*mdl.G/SR;
    d = SQ*(Xref-mdl.O*current_lift(mdl,X)-mdl.G*Uref);
    [left,singular,~] = svd(B,'econ');
    projected = left'*d;
    nullResidual = d-left*projected;
    profile = struct('lambda',diag(singular).^2,'energy',sum(projected.^2,2), ...
        'nullEnergy',sum(nullResidual.^2,'all'),'totalEnergy',sum(d.^2,'all'), ...
        'M',size(X,2),'dimension',size(d,1)*size(d,2),'sensitivity',norm(B,'fro'));
end

function [objective, S, logdet] = profile_value(logRatio, profile)
%profile_value Equation (24), omitting constants independent of the ratio.
    rho = exp(logRatio); terms = rho*profile.lambda;
    S = profile.nullEnergy+sum(profile.energy./(1+terms));
    logdet = sum(log1p(terms));
    if ~isfinite(S) || S <= 0, objective = -Inf; return; end
    objective = -profile.dimension/2*log(S)-profile.M/2*logdet;
end

function [q, r, diagnostic] = maximize_profile(profile, cfg, oldQ, oldR)
%maximize_profile Log grid, every sampled local peak, and explicit end points.
% This is a finite-interval numerical search, not a global-optimum theorem.
    q = oldQ; r = oldR;
    diagnostic = struct('valid',false,'boundary',0,'S',NaN,'objective',NaN,'gain',NaN);
    if profile.totalEnergy <= cfg.energyFloor || profile.sensitivity <= cfg.onlineSensitivityFloor, return; end
    bounds = log(cfg.ratioBounds);
    grid = linspace(bounds(1),bounds(2),cfg.ratioGridSize);
    values = arrayfun(@(x) profile_value(x,profile),grid);
    [best,bestIndex] = max(values); theta = grid(bestIndex);
    opt = optimset('Display','off','TolX',cfg.ratioLogTolerance,'MaxFunEvals',100);
    peaks = find(values(2:end-1) >= values(1:end-2) & values(2:end-1) >= values(3:end))+1;
    for j = peaks
        [trial,negative] = fminbnd(@(x) -profile_value(x,profile),grid(j-1),grid(j+1),opt);
        if -negative > best, theta = trial; best = -negative; end
    end
    previousTheta = min(max(log(oldQ/oldR),bounds(1)),bounds(2));
    previousValue = profile_value(previousTheta,profile);
    if previousValue > best, theta = previousTheta; best = previousValue; end
    [~,S] = profile_value(theta,profile);
    proposedQ = profile.dimension/(2*S); proposedR = proposedQ/exp(theta);
    if ~all(isfinite([proposedQ,proposedR])) || min(proposedQ,proposedR) <= 0, return; end
    q = proposedQ; r = proposedR; diagnostic.valid = true;
    diagnostic.S = S; diagnostic.objective = best; diagnostic.gain = best-previousValue;
    if abs(theta-bounds(1)) <= 10*cfg.ratioLogTolerance, diagnostic.boundary = -1; end
    if abs(theta-bounds(2)) <= 10*cfg.ratioLogTolerance, diagnostic.boundary = 1; end
end

function out = simulate_mpc(mdl, cfg, weights, online, x0, xref)
%simulate_mpc Algorithm 2 updates model and weights BEFORE the current MPC move.
    N = cfg.N; q = weights(1); r = weights(2);
    out = struct('x',zeros(2,N+1),'u',zeros(1,N),'q',zeros(N,1),'r',zeros(N,1), ...
        'uref',evaluation_input_reference(xref(:,1:N),cfg), ...
        'controllerUref',zeros(1,N),'held',false(N,1),'paused',false(N,1), ...
        'trackingRMS',nan(N,1),'resumeCount',0,'stepSeconds',zeros(N,1), ...
        'modelAccepted',false(N,1),'modelRejected',false(N,1),'modelFitError',nan(N,1), ...
        'ratioBoundary',zeros(N,1),'profileS',nan(N,1),'profileObjective',nan(N,1), ...
        'profileGain',nan(N,1),'referenceDegenerate',false(N,1),'qpFallback',false(N,1));
    out.x(:,1) = x0;
    if online
        [controlModel,identification] = initialize_online_model(mdl,cfg);
        out.initialW = identification.W0;
        nZ = size(controlModel.A,1); out.modelA = zeros(nZ,nZ,N); out.modelB = zeros(nZ,N);
    else
        controlModel = mdl;
        [fixedInputReference,~] = make_input_reference(mdl,xref,cfg);
    end
    windowX = zeros(2,cfg.window); windowRef = zeros(2*cfg.Np,cfg.window);
    windowInputRef = zeros(cfg.Np,cfg.window); windowError = zeros(2,cfg.window);
    count = 0; paused = false; previousInput = 0;
    for k = 1:N
        timer = tic;
        if online && k >= 2
            [controlModel,identification,accepted,fitError] = recursive_online_update( ...
                controlModel,identification,out.x(:,k-1),out.u(k-1),out.x(:,k),cfg);
            out.modelAccepted(k) = accepted; out.modelRejected(k) = ~accepted;
            out.modelFitError(k) = fitError;
        end
        if online
            [inputReference,referenceInfo] = make_input_reference(controlModel,xref(:,k),cfg);
            out.referenceDegenerate(k) = referenceInfo.degenerate;
        else
            inputReference = fixedInputReference(k);
        end
        future = repmat(xref(:,k),cfg.Np,1); futureInput = repmat(inputReference,cfg.Np,1);
        if online
            if count < cfg.window
                count = count+1;
            else
                windowX(:,1:end-1) = windowX(:,2:end); windowRef(:,1:end-1) = windowRef(:,2:end);
                windowInputRef(:,1:end-1) = windowInputRef(:,2:end); windowError(:,1:end-1) = windowError(:,2:end);
            end
            windowX(:,count) = out.x(:,k); windowRef(:,count) = future;
            windowInputRef(:,count) = futureInput; windowError(:,count) = out.x(:,k)-xref(:,k);
            if count == cfg.window
                rmsError = sqrt(mean(sum(windowError.^2,1))); out.trackingRMS(k) = rmsError;
                if cfg.enableOnlinePause
                    if paused && max(rmsError,norm(windowError(:,count))) >= cfg.onlineResumeError
                        paused = false; out.resumeCount = out.resumeCount+1;
                    elseif ~paused && rmsError <= cfg.onlinePauseError
                        paused = true;
                    end
                end
                if ~paused
                    profile = prepare_profile(controlModel,windowX,windowRef,windowInputRef);
                    [q,r,diagnostic] = maximize_profile(profile,cfg,q,r);
                    out.held(k) = ~diagnostic.valid;
                    out.ratioBoundary(k) = diagnostic.boundary; out.profileS(k) = diagnostic.S;
                    out.profileObjective(k) = diagnostic.objective; out.profileGain(k) = diagnostic.gain;
                else
                    out.held(k) = true;
                end
            end
            out.paused(k) = paused; out.modelA(:,:,k) = controlModel.A; out.modelB(:,k) = controlModel.B;
        end
        try
            U = mpc_batch(controlModel,out.x(:,k),future,q,r,futureInput);
            controlInput = U(1);
        catch exception
            if ~online, rethrow(exception); end
            controlInput = min(max(previousInput,-cfg.uMax),cfg.uMax); out.qpFallback(k) = true;
        end
        out.u(k) = controlInput; out.q(k) = q; out.r(k) = r;
        out.controllerUref(k) = inputReference; out.stepSeconds(k) = toc(timer);
        out.x(:,k+1) = current_plant_step(mdl,out.x(:,k),controlInput,cfg,k);
        assert(all(isfinite(out.x(:,k+1))), 'Nonfinite true plant state.');
        previousInput = controlInput;
    end
    out.Jx = sum(vecnorm(out.x(:,1:N)-xref(:,1:N),2,1));
    out.Ju = sum(abs(out.u-out.uref));
    out.terminalError = norm(out.x(:,end)-xref(:,N+1));
    out.heldUpdates = nnz(out.held); out.meanStepSeconds = mean(out.stepSeconds);
    if online
        out.finalIdentification = identification;
        out.transitionCount = identification.samples;
    end
end

function write_current_diagnostics(runs,names,cfg,outputFolder,mdl,offlineX,offlineXref,offlineUref,offlineTrace)
%write_current_diagnostics Save quantities required to review Algorithm 2.
    s = runs{4}; t = (0:cfg.N-1)'*cfg.dt;
    diagnostics = table(t,s.q,s.r,s.q./s.r,s.controllerUref',s.uref',s.paused,s.held, ...
        s.modelAccepted,s.modelRejected,s.modelFitError,s.ratioBoundary,s.profileS, ...
        s.profileObjective,s.profileGain,s.referenceDegenerate,s.qpFallback,s.stepSeconds, ...
        'VariableNames',{'Time','q','r','WeightRatio','ControllerInputReference','EvaluationInputReference', ...
        'Paused','Held','ModelAccepted','ModelRejected','OneStepPredictionError','RatioBoundary', ...
        'ProfileS','WindowProfileObjective','ProfileImprovement','ReferenceDegenerate','QpFallback','StepSeconds'});
    writetable(diagnostics,fullfile(outputFolder,'online_diagnostics.csv'));
    save(fullfile(outputFolder,'online_model_history.mat'),'s','-v7');
    profile = prepare_profile(mdl,offlineX,offlineXref,offlineUref);
    [qp,rp,search] = maximize_profile(profile,cfg,offlineTrace(end,2),offlineTrace(end,3));
    fixedPointRatio = offlineTrace(end,2)/offlineTrace(end,3);
    fixedPointValue = profile_value(log(fixedPointRatio),profile);
    offlineVerification = table(offlineTrace(end,2)/offlineTrace(end,3),qp/rp,search.boundary, ...
        search.objective-fixedPointValue, ...
        fixedPointRatio < cfg.ratioBounds(1) || fixedPointRatio > cfg.ratioBounds(2), ...
        'VariableNames',{'FixedPointRatio','IntervalProfileRatio','ProfileBoundary', ...
        'ProfileGainOverFixedPointRatio','FixedPointOutsideRatioInterval'});
    writetable(offlineVerification,fullfile(outputFolder,'offline_profile_crosscheck.csv'));
    summary = zeros(numel(runs),8);
    for j = 1:numel(runs)
        run = runs{j};
        if isfield(run,'stepSeconds')
            times = sort(run.stepSeconds);
            summary(j,1:3) = [mean(times),times(max(1,ceil(0.95*numel(times)))),max(times)];
            summary(j,4) = nnz(times > cfg.dt);
        else
            summary(j,1:3) = [run.meanStepSeconds,NaN,NaN]; summary(j,4) = NaN;
        end
        summary(j,5) = max([0,abs(run.u)-cfg.uMax]);
        if isfield(run,'ratioBoundary')
            summary(j,6:8) = [nnz(run.ratioBoundary),nnz(run.modelRejected),nnz(run.qpFallback)];
        end
    end
    timing = array2table(summary,'VariableNames',{'MeanSeconds','P95Seconds','MaxSeconds', ...
        'Overruns','MaximumInputViolation','BoundaryUpdates','RejectedModels','QpFallbacks'});
    timing = addvars(timing,names','Before',1,'NewVariableNames','Controller');
    writetable(timing,fullfile(outputFolder,'timing_and_guards.csv'));
    fid = fopen(fullfile(outputFolder,'algorithm_settings.json'),'w','n','UTF-8');
    assert(fid >= 0); cleanup = onCleanup(@() fclose(fid));
    fprintf(fid,'%s\n',jsonencode(cfg,PrettyPrint=true));
    assert(all(isfinite(s.x(:))) && all(isfinite(s.q)) && all(isfinite(s.r)) ...
        && all(s.q > 0) && all(s.r > 0));
    assert(s.transitionCount == cfg.N-1,'A measured transition was omitted.');
    assert(all(s.profileGain(isfinite(s.profileGain)) >= -1e-8),'Profile search decreased its window objective.');
end

function verification = verify_current_algorithm(cfg)
%verify_current_algorithm Check actual public run-path helpers before simulation.
    localCfg = cfg; localCfg.Np = 1; localCfg.Q0 = diag([2,3]); localCfg.R0 = 4;
    localCfg.uMax = 0.25;
    model = struct('A',eye(2),'B',[1;0],'C',eye(2),'O',zeros(2),'G',[1;0], ...
        'Qbar',localCfg.Q0,'Rbar',4,'uMax',0.25,'K',2,'GQ',[2,0], ...
        'qpOptions',optimoptions('quadprog','Display','off','Algorithm','interior-point-convex'));
    X = zeros(2,2); references = [3,1;0.2,-0.4]; Uref = zeros(1,2);
    [actual,stats] = mpc_batch(model,X,references,3,2,Uref,true);
    [virtual,virtualStats] = mpc_batch(model,X,references,3,2,Uref,false);
    profile = prepare_profile(model,X,references,Uref);
    [~,S,ld] = profile_value(log(3/2),profile);
    logdetR = log(4); logdetH = log(2*(3*2+2*4));
    gProfile = profile.dimension/2*log(3)-3*S-profile.M/2*ld;
    expected = gProfile-profile.M/2*(logdetR+log(2));
    verification.evidenceIdentityError = abs(stats.score-expected);
    verification.virtualCostIdentityError = abs(3*stats.EX+2*stats.EU-3*S);
    verification.constraintGap = norm(actual-virtual,'fro');
    verification.virtualStatsIndependentOfConstraints = abs(stats.score-virtualStats.score);
    verification.hessianLogdet = logdetH;
    qStep = 3e-5; rStep = 2e-5;
    [~,qp] = mpc_batch(model,X,references,3+qStep,2,Uref,false);
    [~,qm] = mpc_batch(model,X,references,3-qStep,2,Uref,false);
    [~,rp] = mpc_batch(model,X,references,3,2+rStep,Uref,false);
    [~,rm] = mpc_batch(model,X,references,3,2-rStep,Uref,false);
    analytic = [(profile.M*(2-1)+2*2*stats.h)/(2*3)-stats.EX; profile.M/(2*2)-stats.EU-stats.h];
    numeric = [(qp.score-qm.score)/(2*qStep);(rp.score-rm.score)/(2*rStep)];
    verification.gradientError = max(abs(analytic-numeric));
    [qStar,rStar,search] = maximize_profile(profile,localCfg,3,2);
    verification.positiveProfileWeights = search.valid && min(qStar,rStar) > 0;
    zero = model; zero.B(:) = 0;
    [fallback,info,unchanged] = make_input_reference(zero,[1;0],localCfg);
    verification.zeroSensitivityGuard = isfinite(fallback) && info.degenerate && isequal(unchanged,[1;0]);
    localCfg.onlineA0 = []; localCfg.onlineB0 = [];
    [cold,identification] = initialize_online_model(model,localCfg);
    [updated,identification,accepted] = recursive_online_update(cold,identification, ...
        [1;0.2],0.4,[0.9;0.1],localCfg);
    batchW = identification.accumulator/identification.information;
    verification.recursiveBatchError = norm([updated.A,updated.B]-batchW,'fro');
    unrelated = model; unrelated.A = 10*eye(2); unrelated.B = [-999;999];
    [~,otherIdentification] = initialize_online_model(unrelated,localCfg);
    verification.initialModelIndependentOfOfflineFit = isequal( ...
        otherIdentification.W0,[cold.A,cold.B]);
    assert(verification.evidenceIdentityError < 1e-9 && verification.virtualCostIdentityError < 1e-9);
    assert(verification.gradientError < 1e-7 && verification.constraintGap > 0);
    assert(verification.virtualStatsIndependentOfConstraints < 1e-12);
    assert(verification.positiveProfileWeights && verification.zeroSensitivityGuard);
    assert(accepted && verification.recursiveBatchError < 1e-9);
    assert(verification.initialModelIndependentOfOfflineFit);
end

function next = current_plant_step(~,x,u,cfg,k)
%current_plant_step Simulate the physical plant only after applying the move.
    gain = 1;
    if (k-1)*cfg.dt >= cfg.gainChangeTime, gain = cfg.gainAfterChange; end
    next = quench_step(x,u,cfg,gain);
end

function uref = evaluation_input_reference(xref,cfg)
%evaluation_input_reference Common physical equilibrium for METRICS ONLY.
% This reference is never passed to the online controller or its evidence.
    h = cfg.hNom+cfg.hScale*xref(1,:);
    head0 = cfg.pressure-cfg.downstreamPressure+cfg.rhoWater*cfg.g*cfg.hNom;
    factor = sqrt((cfg.pressure-cfg.downstreamPressure+cfg.rhoWater*cfg.g*h)/head0);
    uref = (cfg.massNom./factor-cfg.massNom)/cfg.massScale;
    changed = (0:size(xref,2)-1)*cfg.dt >= cfg.gainChangeTime;
    uref(changed) = uref(changed)/cfg.gainAfterChange;
end
