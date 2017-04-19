clc
clear all
close all

s = RandStream('mt19937ar', 'Seed', 1);
RandStream.setGlobalStream(s);
addpath(genpath('../../../toolboxes'));
% Generate Samples from the model
N = 100;
t = linspace(0, 5, N)';
D = 2;
kern = kernCreate(t, {'multi', 'rbf', 'sim', 'sim'});
lengthscale = 1;
decay = [0.4 1];
sensitivityc = [1 1];
% Initialize the lengthscale for the rbf kernel
kern.comp{1}.inverseWidth = 2/(lengthscale.^2);
for d=1:D,
   kern.comp{1+d}.inverseWidth = 2/(lengthscale.^2);
   kern.comp{1+d}.decay = decay(d);
   kern.comp{1+d}.sensitivity = sensitivityc(d);
end
K = kernCompute(kern, t);
yV = real(gsamp(zeros(size(K, 1),1), K, 1))';
% Get the data and plot it
Xt = cell(1, D);
yt = cell(1, D);
startVal = N+1;
endVal = N;
figure
for d = 1:D
    endVal = endVal + N;
    yt{d} = yV(startVal:endVal);
    Xt{d} = t;
    startVal = endVal + 1;
    subplot(2,1,d)
    plot(Xt{d}, yt{d},'linewidth',1.5)
    xlabel('Time')
    ylabel(['Output ' num2str(d)])
end
% Get the training set and the test set
X = cell(1, D + 1);
y = cell(1, D + 1);
XTest = cell(1, D + 1);
yTest = cell(1, D + 1);
% Necessary for the structure of the kernel in MATLAB
options = multigpOptions('ftc');
options.kernType = 'sim';
options.optimiser = 'scg';
options.nlf = 1;

X{1} = 0;
XTest{1} = 0;
yTrain{1} = [];
nTrain = 50;
for d=1:D
    index = randperm(N);
    indexTrain = sort(index(1:nTrain));
    indexTest = sort(index(nTrain+1:end));
    X{d+1} = Xt{d}(indexTrain);
    y{d+1} = yt{d}(indexTrain);
%     options.bias(d+1) = mean(y{d+1});
%     options.scale(d+1) = std(y{d+1});
    XTest{d+1} = Xt{d}(indexTest);
    yTest{d+1} = yt{d}(indexTest);
end
figure
for d = 1:D
    endVal = endVal + N;    
    startVal = endVal + 1;
    subplot(2,1,d)
    hold on
    plot(Xt{d}, yt{d}, 'linewidth', 1.5)
    hold on
    plot(X{d+1}, y{d+1}, '.r', 'markersize', 20)
    xlabel('Time')
    ylabel(['Output ' num2str(d)])
    legend('All data', 'Training data')
end

display = 2;
iters = 200;
q = size(X{1}, 2);
d = D + options.nlf;
model = multigpCreate(q, d, X, y, options);
% Change initial parameters
[paramst,names] = modelExtractParam(model);
% index = paramNameRegularExpressionLookup(model, '.* white .* variance');
% params(index) = log(10);
%index = paramNameRegularExpressionLookup(model, '.* inverse width latent .*');
%params(index) = [log(1) log(10)];
%model = modelExpandParam(model, params);
%Trains the model
model = multigpOptimise(model, display, iters);
% Compute the posterior mean
[mu, var] = multigpPosteriorMeanVar(model, Xt{1});


% Plot the predictions,
figure
for d = 1:D,
    endVal = endVal + N;
    startVal = endVal + 1;
    subplot(2,1,d)
    hold on
    plot(Xt{d}, yt{d}, 'linewidth', 1.5)
    hold on
    plot(X{d+1}, y{d+1}, '.r', 'markersize', 20)
    hold on
    plot(Xt{d}, mu{d+1}, 'k', 'linewidth', 1.5)
    xlabel('Time')
    ylabel(['Output ' num2str(d)])
    legend('All data', 'Training data', 'Mean prediction')
end
[params, names] = modelExtractParam(model);