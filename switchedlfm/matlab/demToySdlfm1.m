% DEMTOYSDLFM1 Demonstrate switching dynamical latent force model on
% artificial data. Reproduces one try of toy example 1 in the paper Álvarez
% et al (2011): Switched Latent Force Models for Movement Segmentation. at
% NIPS 2010

% SDLFMGP

clc
clear
randn('state', 1e6)
rand('twister', 1e6)

addToolboxes(0,1)

dataSetName = 'toySdlfm1';
experimentNo = 1;

[XTemp2, yTemp2, XTestTemp, yTestTemp] = mapLoadData(dataSetName);

options = sdlfmgpOptions('ftc');
switchingTimes = [-1 5 8];

iters = 200;
display = 1;

options.nIntervals = 3;
options.nlfPerInt = 1;
options.kern.nIntervals = options.nIntervals;
options.kern.switchingTimes = switchingTimes;
options.kern.nlfPerInt = options.nlfPerInt;
options.kern.isNormalised = true;

yTemp = yTemp2{1};
XTemp = XTemp2{1};

Y = cell(1, length(yTemp));
X = cell(1, length(yTemp));
scaleVal = zeros(1, length(yTemp));

for i=1:length(yTemp)
    Y{i} = yTemp{i};
    X{i} = XTemp{i};
    options.scale(i) = std(Y{i});
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
count = length(model.fix);
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
% Change values of initial parameters

% Change values of initial parameters
params = modelExtractParam(model);
indexDamper = paramNameRegularExpressionLookup(model, '.* damper');
params(indexDamper) = log(0.4);
indexSpring = paramNameRegularExpressionLookup(model, '.* spring');
params(indexSpring) = log(5);
indexVariance = paramNameRegularExpressionLookup(model, '.* variance');
params(indexVariance) = log(50);
model = modelExpandParam(model, params);

% Train the model

model = modelOptimise(model, [], [], display, iters);

% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');

[~, ~, XGT, fGT] = mapLoadData(dataSetName);
sdlfmgpToyResults(dataSetName, experimentNo, XTemp, yTemp, ...
    XGT{1}, fGT{1})


