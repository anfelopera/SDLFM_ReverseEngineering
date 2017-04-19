% DEMTOYSDSIM3 Demonstrate switching dynamical single input motif on
% artificial data. Based on the paper Álvarez et al (2011): Switched 
% Latent Force Models for Movement Segmentation. at NIPS 2010

% SDLFMGP

% clc
clear
randn('state', 1e6)
rand('twister', 1e6)

% % add toolboxes
addToolboxes(0,1)

% % defines the model
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

for i=1:length(yTemp)
    Y{i} = yTemp{i};
    X{i} = XTemp{i};
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

capName = dataSetName;
capName(1) = upper(capName(1));
try
    % % load model if the experiment was previously saved 
    load(['dem' capName num2str(experimentNo) '.mat']);
catch
    % Change values of initial parameters
    [params,names]        = modelExtractParam(model);
    indexSensitivity      = paramNameRegularExpressionLookup(model, '.* sensitivity .*');
    params(indexSensitivity) = ([1 1 1 1, 1 1 1 1, 1 1 1 0.1]);
    indexPriorSvariance   = paramNameRegularExpressionLookup(model, '.* priorSvariance .*');
    params(indexPriorSvariance) = log(0.3);
    indexVariance         = paramNameRegularExpressionLookup(model, '.* variance');
    params(indexVariance) = log(1e-3);
    model = modelExpandParam(model, params);

    % Train the model
    model = modelOptimise(model, [], [], display, iters);

    % Save the results.
    save(['dem' capName num2str(experimentNo) '.mat'], 'model');
end

% % results
[~, ~, XGT, fGT] = mapLoadData(dataSetName);
sparseSdsimgpToyResults(dataSetName, experimentNo, XTemp, yTemp, XGT{1}, fGT{1})
[param,names] = modelExtractParam(model);

indexSP = paramNameRegularExpressionLookup(model, '.* switching point .*');
for k =1:d+1
    figure(k); hold on; 
    estimatedSwitchingTimes  = exp(param(indexSP(2:end)));
    estimatedSwitchingTimes  = [switchingTimes(1) estimatedSwitchingTimes];
    lines = cumsum(estimatedSwitchingTimes);
    for j = 1:options.nIntervals-1
        plot(lines(j+1)*[1 1],[-50 50],'k--','linewidth',1)
    end
end
