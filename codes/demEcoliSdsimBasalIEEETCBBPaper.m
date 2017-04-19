% demEcoliSdsimBasalIEEETCBB Demonstrate switching dynamical single input
% motif on Ecoli data from:
%   Sanguinetti et al., "Switching regulatory models of cellular stress
%   response", 2009.
% Based on the paper Álvarez et al (2011): Switched Latent Force Models
% for Movement Segmentation. at NIPS 2010

% clc
clear
randn('state', 1e6)
rand('twister', 1e6)

addToolboxes(0,1)

dataSetName = 'toyEcoliSdsim';
experimentNo = 1;

load('dataEcoliMicroData.mat')
genes2 = {'hypB','aspA','yciD','yjiD','moaA'};

for g = 1:length(genes2)
    gene{g} = genes2{g};
    idx = find(strcmp(gene{g},genes));        
    XTemp2{1}{g} = times';
    yTemp2{1}{g} = [1;data(idx,:)'];
    XTestTemp{1}{g} = (0:0.1:times(end)+1)';
    XTestTemp{1}{g}(times(2:end) +1 ) = [];
    yTestTemp{1}{g} = [];
end

options                     = sdlfmgpOptions('ftc');
options.nIntervals          = 2;
switchingTimes              = 16.6;%21.4;
if options.nIntervals == 1
    switchingTimes          = -1;
else
    switchingTimes          = [-1 switchingTimes ];
end

iters                       = 200;
display                     = 1;

options.nlfPerInt           = 1;
options.kern.nIntervals     = options.nIntervals;
options.kern.switchingTimes = switchingTimes;
options.kern.nlfPerInt      = options.nlfPerInt;
options.kern.isNegativeS    = true;
options.kern.isNormalised   = true;
options.kern.isStationary   = false;
options.kern.priorS         = false;
options.kern.delay          = 0;
options.kernType            = 'sdsim'; 

yTemp = yTemp2{1};
XTemp = XTemp2{1};

Y = cell(1, length(yTemp));
X = cell(1, length(yTemp));
scaleVal = zeros(1, length(yTemp));

for i=1:length(yTemp)
    Y{i} = yTemp{i};
    X{i} = XTemp{i};
%     options.scale(i) = std(Y{i});
end

q = 1;
d = length(yTemp);

% Creates the model
ttime = cputime;
options.X = X;
options.meanFunction = meanCreate(q, d, options);
model = sdlfmgpCreate(q, d, X, Y, options);

% % We fix some additional parameters, namely, the value of the first
% % switching point at -1 and the value of the covariances for the initial
% % conditions.
model2 = model;
model2.nParams = size(model2.paramGroups,1);
model2.paramGroups = speye(model2.nParams);
count = 0;
% % Fix the value of the initial condition first
index = paramNameRegularExpressionLookup(model2, '.* switching point interval 1');
for k=1:length(index);
    count = count + 1;
    model.fix(count).index = index(k);
    model.fix(count).value = model.kern.comp{1}.comp{1}.switchingTimes(1);
end
% % Fix the value of the second condition
index = paramNameRegularExpressionLookup(model2, '.* switching point interval 2');
for k=1:length(index);
    count = count + 1;
    model.fix(count).index = index(k);
    model.fix(count).value = log(model.kern.comp{1}.comp{1}.switchingTimes(2));
end

% % Now, fix the covariances of the initial conditions
valkInit = model.kern.comp{1}.LIC(:)';
index = paramNameRegularExpressionLookup(model2, '.* kInit.*');
for k=1:length(index);
    count = count + 1;
    model.fix(count).index = index(k);
    model.fix(count).value = valkInit(k);
end

% % Tie the basals and decays parameters
indexDecay = paramNameRegularExpressionLookup(model, '.* decay');
indexBasal = paramNameRegularExpressionLookup(model, '.* basal');
indexBasal2 = paramNameRegularExpressionLookup(model2, '.* basal');
for i = d:-1:1
    model.paramGroups(indexBasal2(i),indexDecay(i)) = 1;
    model.paramGroups(:,indexBasal(i)) = [];
end
% % Force the sesitivity to be positive
for i = 1:d
    model.kern.comp{1,1}.comp{1,i}.transforms.index = ...
        [model.kern.comp{1,1}.comp{1,i}.transforms.index, ... 
        model.kern.comp{1,1}.comp{1,i}.sensitivityIndx];
end
% % Tie the sensitivity per each gene
for i = d:-1:1
    indexSensitivity = paramNameRegularExpressionLookup(model, [ '.*' num2str(i) ' sensitivity .*']);
    indexSensitivity2 = paramNameRegularExpressionLookup(model2, [ '.*' num2str(i) ' sensitivity .*']);
    for k = length(switchingTimes):-1:2          
        model.paramGroups(indexSensitivity2(k),indexSensitivity(1)) = 1;
        model.paramGroups(:,indexSensitivity(k)) = [];         
    end
end
% % % Tie the noise variance
% indexVariance = paramNameRegularExpressionLookup(model, '.* variance');
% indexVariance2 = paramNameRegularExpressionLookup(model2, '.* variance');
% for i = d:-1:2
%     model.paramGroups(indexVariance2(i),indexVariance(1)) = 1;
%     model.paramGroups(:,indexVariance(i)) = [];
% end

capName = dataSetName;
capName(1) = upper(capName(1));
try
    % % load model if the experiment was previously saved 
    load(['dem' capName num2str(experimentNo) '.mat']);
catch
    % % Change values of initial parameters
    [params,names] = modelExtractParam(model);
    indexDecay = paramNameRegularExpressionLookup(model, '.* decay');
    params(indexDecay) = log([0.15 0.15 0.05 0.03 0.15]);
    indexInverseW = paramNameRegularExpressionLookup(model, '.* inverse width');
    params(indexInverseW) = log([1e-2 1e-4]);
    indexSensitivity = paramNameRegularExpressionLookup(model, '.* sensitivity .*');
    params(indexSensitivity) = log([0.4 0.2 1.3 0.25 0.30]);
    indexVariance = paramNameRegularExpressionLookup(model, '.* variance');
    params(indexVariance) = log([1e-2 0.5e-2 1e-1 2e-2 1e-2]);
    model = modelExpandParam(model, params);

    % Train the model
    model = modelOptimise(model, [], [], display, iters);

    % Save the results.
    save(['dem' capName num2str(experimentNo) '.mat'], 'model');
end

XGT = XTestTemp;
fGT = yTestTemp;
sdsimgpToyResultsEcoli(dataSetName, experimentNo, X, Y, ...
    XGT{1}, fGT{1})
load(['dem' capName num2str(experimentNo) '.mat']);
[param,names] = modelExtractParam(model);
indexSP = paramNameRegularExpressionLookup(model, '.* switching point .*');

for k =1:d+options.nlfPerInt
    figure(k); hold on; 
    estimatedSwitchingTimes  = [-1 exp(param(indexSP(2:end)))];
    lines = cumsum(estimatedSwitchingTimes);
    for j = 1:options.nIntervals-1
        plot(lines(j+1)*[1 1],[-50 50],'k--','linewidth',1)
    end
end

figure(1); ylim([0 1])
figure(2); ylim([0.5 4.5])
figure(3); ylim([0.8 2.8])
figure(4); ylim([0 35])
figure(5); ylim([0 3])
figure(6); ylim([0.5 3.5])
ttime = cputime-ttime
