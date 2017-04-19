% DEMTOYSDSIM1 Demonstrate switching dynamical single input motif on
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
dataSetName  = 'toySdsim1';
experimentNo = 1;

[XTemp2, yTemp2, XTestTemp, yTestTemp] = mapLoadData(dataSetName);

options                     = sdlfmgpOptions('ftc');
switchingTimes              = [-1 6 7];

iters                       = 50;
display                     = 1;

options.nIntervals          = 3;
options.nlfPerInt           = 1;
options.kern.nIntervals     = options.nIntervals;
options.kern.switchingTimes = switchingTimes;
options.kern.nlfPerInt      = options.nlfPerInt;
options.kern.isNegativeS    = true;
options.kern.isNormalised   = true;
options.kern.isStationary   = false;
options.kern.delay          = 0;
options.kernType            = 'sdsim'; 

yTemp = yTemp2{1};
XTemp = XTemp2{1};

Y = cell(1, 1);
X = cell(1, 1);
Y{1} = yTemp;
X{1} = XTemp;

% Get the time index.
q = 1;
d = 1;

% Creates the model
model = sdlfmgpCreate(q, d, X, Y, options);

% % We fix some additional parameters, namely, the value of the first
% % switching point at -1 and the value of the covariances for the initial
% % conditions.
model2 = model;
model2.nParams = size(model2.paramGroups,1);
model2.paramGroups = speye(model2.nParams);
count = 0;
% % Fix the value of the initial condition first
% % Now, fix the covariances of the initial conditions
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

capName = dataSetName;
capName(1) = upper(capName(1));
try
    % % load model if the experiment was previously saved 
    load(['dem' capName num2str(experimentNo) '.mat']);
catch
    % % Change values of initial parameters
    [params,names]= modelExtractParam(model);
    indexDecay            = paramNameRegularExpressionLookup(model, '.* decay');
    params(indexDecay)    = log(1);
    indexVariance         = paramNameRegularExpressionLookup(model, '.* variance');
    params(indexVariance) = log(0.01);
    model = modelExpandParam(model, params);

    % Train the model
    model = modelOptimise(model, [], [], display, iters);
    
    % Save the results.
    save(['dem' capName num2str(experimentNo) '.mat'], 'model');
end

% % results
[~, ~, XGT, fGT] = mapLoadData(dataSetName);
sdsimgpToyResults(dataSetName, experimentNo, X, Y, XGT{1}, fGT{1})

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
