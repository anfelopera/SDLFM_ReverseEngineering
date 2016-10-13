% DEMTOYSDSIM3 Demonstrate switching dynamical single input motif on
% artificial data. Based on the paper Álvarez et al (2011): Switched 
% Latent Force Models for Movement Segmentation. at NIPS 2010

% SDLFMGP

% clc
clear
randn('state', 1e6)
rand('twister', 1e6)

addToolboxes(0,1)

dataSetName  = 'toySparseSdsim1';
experimentNo = 2;

[XTemp2, yTemp2, XTestTemp, yTestTemp] = mapLoadData(dataSetName);

options                     = sdlfmgpOptions('ftc');
switchingTimes              = [-1 10];

iters                       = 200;
display                     = 1;

options.nIntervals          = 2;
options.nlfPerInt           = 2;
options.kern.nIntervals     = options.nIntervals;
options.kern.switchingTimes = switchingTimes;
options.kern.nlfPerInt      = options.nlfPerInt;
options.kern.isNegativeS    = true;
options.kern.isNormalised   = true;
options.kern.isStationary   = false;
options.kern.priorS         = true;
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

% Get the time index.

q = 1;
d = length(yTemp);

% Creates the model

model = sdlfmgpCreate(q, d, X, Y, options);

% We fix some additional parameters, namely, the value of the first
% switching point at -1 and the value of the covariances for the initial
% conditions.
model2 = model;
model2.nParams = size(model2.paramGroups,1);
model2.paramGroups = speye(model2.nParams);
count = 0;
% Fix the value of the initial condition first
% Now, fix the covariances of the initial conditions
valkInit = model.kern.comp{1}.LIC(:)';
index = paramNameRegularExpressionLookup(model2, '.* kInit.*');
for k=1:length(index);
    count = count + 1;
    model.fix(count).index = index(k);
    model.fix(count).value = valkInit(k);
end
index = paramNameRegularExpressionLookup(model2, '.* switching point interval 1');
for k=1:length(index);
    count = count + 1;
    model.fix(count).index = index(k);
    model.fix(count).value = model.kern.comp{1}.comp{1}.switchingTimes(1);
end

index = paramNameRegularExpressionLookup(model2, '.* switching point interval 2');
for k=1:length(index);
    count = count + 1;
    model.fix(count).index = index(k);
    model.fix(count).value = log(model.kern.comp{1}.comp{1}.switchingTimes(2));
end

decayVal = log([0.7 1.5 0.2]);
index = paramNameRegularExpressionLookup(model2, '.* decay');
for k=1:length(index);
    count = count + 1;
    model.fix(count).index = index(k);
    model.fix(count).value = decayVal(k);
end

indexInverseWVal = [1e-1 1 1e-2 1, 1e-1 1 1e-2 1, 1e-1 1 1e-2 1];
index = paramNameRegularExpressionLookup(model2, '.* inverse width');
for k=1:length(index);
    count = count + 1;
    model.fix(count).index = index(k);
    model.fix(count).value = indexInverseWVal(k);
end

% Change values of initial parameters
[params,names]        = modelExtractParam(model);
valkInit              = model.kern.comp{1}.LIC(:)';
indexvalkInit         = paramNameRegularExpressionLookup(model, '.* kInit.*');
params(indexvalkInit) = valkInit;
indexSP1              = paramNameRegularExpressionLookup(model, '.* switching point interval 1');
params(indexSP1)      = model.kern.comp{1}.comp{1}.switchingTimes(1);
indexSP2              = paramNameRegularExpressionLookup(model, '.* switching point interval 2');
params(indexSP2)      = log(model.kern.comp{1}.comp{1}.switchingTimes(2));
indexDecay            = paramNameRegularExpressionLookup(model, '.* decay');
% params(indexDecay)    = log([0.7 1.5 0.2]);
indexInverseW         = paramNameRegularExpressionLookup(model, '.* inverse width');
% params(indexInverseW) = log([1e-1 1; 1e-2 1]);
indexSensitivity      = paramNameRegularExpressionLookup(model, '.* sensitivity .*');
params(indexSensitivity) = ([1 1 1 1, 1 1 1 1, 1 1 1 0.1]);
indexPriorSvariance   = paramNameRegularExpressionLookup(model, '.* priorSvariance .*');
params(indexPriorSvariance) = log(0.3);
indexVariance         = paramNameRegularExpressionLookup(model, '.* variance');
params(indexVariance) = log([0.0146  0.0058 0.0946]);
params(indexVariance) = log(1e-3);

model = modelExpandParam(model, params);

% Train the model
model = modelOptimise(model, [], [], display, iters);
[param] = modelExtractParam(model);

% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');

[~, ~, XGT, fGT] = mapLoadData(dataSetName);
sparseSdsimgpToyResults(dataSetName, experimentNo, XTemp, yTemp, ...
    XGT{1}, fGT{1})

for k =1:d+2
    figure(k); hold on; 
    estimatedSwitchingTimes  = exp(param(options.nlfPerInt*options.nIntervals+3:options.nlfPerInt*options.nIntervals+1+options.nIntervals));
    lines = cumsum([param(options.nlfPerInt*options.nIntervals+2) estimatedSwitchingTimes]);
    for j = 1:options.nIntervals-1
        plot(lines(j+1)*[1 1],[-50 50],'k--','linewidth',1)
    end
end
