% demYeastSdsimIEEETCBB Demonstrate switching dynamical single input
% motif on Yeast data from:
%   Opper et al., "Learning combinatorial transcriptional dynamics from gene
%   expression data", 2010.
% Based on the paper Álvarez et al (2011): Switched Latent Force Models
% for Movement Segmentation. at NIPS 2010

% clc
clear
randn('state', 1e6)
rand('twister', 1e6)

addToolboxes(0,1)

dataSetName  = 'toyYeastSdsim';
experimentNo = 1;

load('dataYeastTuLocal.mat')

% genes2 = {'YLR030W','TKL2'};
% genes2 = {'YLR183C','YLR030W','TKL2','RPL17B','RPS16B','RPL13A'};
% genes2 = {'YLR183C','YLR030W','TKL2'};
% genes2 = {'YLR183C','YLR030W','TKL2','YOR359W','PFK27'};
genes2 = {'YLR183C','YLR030W','TKL2','YOR359W','PFK27','RPL17B','RPS16B','RPL13A'};


for g = 1:length(genes2)
    gene{g} = genes2{g};
    idx     = find(strcmp(gene{g},longNames));    
%     XTemp2{1}{g} = (0:35)';
    XTemp2{1}{g} = (0:20:700)';
    yTemp2{1}{g} = (1.2.^data(idx(1),:))';
    InitVal(g)   = yTemp2{1}{g}(1);
    FinalVal(g)  = yTemp2{1}{g}(end);    
    XTestTemp{1}{g} = (0:700)';  
    for j = 1:length(XTestTemp{1}{g})
        temp = sum(XTestTemp{1}{g}(j) ==  XTemp2{1}{g});
        if temp == 1
            XTestTemp{1}{g}(j) = NaN;
        end
    end
    XTestTemp{1}{g}(isnan(XTestTemp{1}{g})) = [];
    yTestTemp{1}{g} = [];
end

options                     = sdlfmgpOptions('ftc');

options.nIntervals          = 3;
switchingTimes              = [270.3 220.1];
% switchingTimes              = (max(XTestTemp{1}{1})/options.nIntervals )*ones(1,options.nIntervals-1);
if options.nIntervals == 1
    switchingTimes          = [-1];
else
    switchingTimes          = [-1 switchingTimes ];
end
% switchingTimes              = [-1 181 175 175];
% switchingTimes              = [-1 181 220];

iters                       = 20;
display                     = 1;
options.nlfPerInt           = 2;
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

Y = cell(1, length(yTemp));
X = cell(1, length(yTemp));
scaleVal = zeros(1, length(yTemp));

Y{1} = yTemp2{1};
X{1} = XTemp2{1};
% options.scale(1) = std(Y{1});

for i=1:length(yTemp)
    Y{i} = yTemp{i};
    X{i} = XTemp{i};
%     options.scale(i) = std(Y{i});
end

% Get the time index.

q = 1;
d = length(yTemp);

% Creates the model

ttime = cputime;

options.X = X;
options.meanFunction = meanCreate(q, d, options);
% options.priorS = priorSCreate(q, d, options);
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

% % % Fix the value of the second condition first
% index = paramNameRegularExpressionLookup(model2, '.* switching point interval 2');
% for k=1:length(index);
%     count = count + 1;
%     model.fix(count).index = index(k);
%     model.fix(count).value = log(model.kern.comp{1}.comp{1}.switchingTimes(2));
% end
% % % Fix the value of the third condition first
% index = paramNameRegularExpressionLookup(model2, '.* switching point interval 3');
% for k=1:length(index);
%     count = count + 1;
%     model.fix(count).index = index(k);
%     model.fix(count).value = log(model.kern.comp{1}.comp{1}.switchingTimes(3));
% end

% % % Tie the sensitivity per each gene
for i = d:-1:1
    for j = 2:-1:1
    indexSensitivity = paramNameRegularExpressionLookup(model, [ '.*' num2str(i) ' sensitivity ' num2str(j) '.*']);
    indexSensitivity2 = paramNameRegularExpressionLookup(model2, [ '.*' num2str(i) ' sensitivity ' num2str(j) '.*']);
        for k = length(switchingTimes):-1:2
            model.paramGroups(indexSensitivity2(k),indexSensitivity(1)) = 1;
            model.paramGroups(:,indexSensitivity(k)) = [];
        end
    end
end

% Fix the value of the inverse width per each latent force
for i = options.nlfPerInt:-1:1
  indexInverseW = paramNameRegularExpressionLookup(model, ['.* inverse width ' num2str(i) '.*']);
  indexInverseW2 = paramNameRegularExpressionLookup(model2, ['.* inverse width ' num2str(i) '.*']);
  for j = length(indexInverseW):-1:2
    for k = length(indexInverseW2):-1:2
        model.paramGroups(indexInverseW2(k),indexInverseW(1)) = 1;    
    end
    model.paramGroups(:,indexInverseW(j)) = [];
  end
end

% % Tie the noise variance
indexVariance = paramNameRegularExpressionLookup(model, '.* variance');
indexVariance2 = paramNameRegularExpressionLookup(model2, '.* variance');
for i = d:-1:2
    model.paramGroups(indexVariance2(i),indexVariance(1)) = 1;
    model.paramGroups(:,indexVariance(i)) = [];
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
    params(indexDecay)    = log([0.3 1 0.8 0.5 0.5 0.4 0.4 0.4]);
    indexBasal            = paramNameRegularExpressionLookup(model, '.* basal');
    params(indexBasal)    = log([1.2 0.3 1 1.1 0.6 1.2 1 1.2]);
    indexInverseW         = paramNameRegularExpressionLookup(model, '.* inverse width');
    params(indexInverseW) = log([5e-4 7e-4]);
    indexSensitivity      = paramNameRegularExpressionLookup(model, '.* sensitivity .*');
    params(indexSensitivity) = ([-0.4 0,...
                                 0.7 0,...
                                 2.5 0,...
                                 0 1.8,...
                                 0 1.8,...
                                 0.4 2.2,...
                                 0.3 2.1,...
                                 0.3 2.2]);                        
    indexVariance         = paramNameRegularExpressionLookup(model, '.* variance');
    params(indexVariance) = log([1e-2]);
    model = modelExpandParam(model, params);

    % Train the model
    model = modelOptimise(model, [], [], display, iters);

    % Save the results.
    save(['dem' capName num2str(experimentNo) '.mat'], 'model');
end

XGT = XTestTemp;
fGT = yTestTemp;
% [~, ~, XGT, fGT] = mapLoadData(dataSetName);
sdsimgpToyResultsYeast(dataSetName, experimentNo, X, Y, ...
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

figure(1); ylim([0 1]); title('FHL1')
figure(2); ylim([0 1]); title('RAP1')
figure(3); ylim([1 5]); title(gene{1})
figure(4); ylim([0 2.5]); title(gene{2})
figure(5); ylim([0 8]); title(gene{3})
figure(6); ylim([1.5 6]); title(gene{4})
figure(7); ylim([1 5]); title(gene{5})
figure(8); ylim([3 10]); title(gene{6})
figure(9); ylim([2 9]); title(gene{7})
figure(10); ylim([4 10]); title(gene{8})
ttime = cputime-ttime