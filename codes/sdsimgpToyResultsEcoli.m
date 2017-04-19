function [mu,varsigma] = sdsimgpToyResultsEcoli(dataSetName, experimentNo, XTemp, yTemp, XGT, fGT)

% SDSIMGPTOYRESULTS Show the prediction results for the demSdsimgpToy demo.
% FORMAT 
% DESC Show the prediction results for the demSdsimgpToy demo.
% ARG dataSetName : name of the dataset used for the demo
% ARG experimentNo : number of the experiment
% ARG XTemp : input locations training data
% ARG yTemp : output values for training data
%
% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio A. Alvarez, 2015

% SDLFMGP


capName = dataSetName;
capName(1) = upper(capName(1));
load(['dem' capName num2str(experimentNo) '.mat'], 'model');

saveFigures = false;
scaleVal = 1;

fontsize = 25;
linewidth = 2;
markersize = 25;

X = cell(size(yTemp, 2)+model.nlfPerInt,1);
y = cell(size(yTemp, 2)+model.nlfPerInt,1);

for j=1:model.nlfPerInt
   y{j} = 0;
   X{j} = 0;
end
for i = 1:size(yTemp, 2)
  y{i+model.nlfPerInt} = yTemp{i};
  X{i+model.nlfPerInt} = XTemp{i};
end

Xt = XGT{1};
[mu, varsigma] = sdsimgpPosteriorMeanVar(model, Xt);

xlim = [Xt(1) Xt(end)];

close all
nFigs = model.nout+model.nlfPerInt;
mu{1} = (mu{1}-min(mu{1}))/(max(mu{1})-min(mu{1}));

for k=1:nFigs,
    figure
    hold on
%     mu{k} = smooth(mu{k},0.1,'lowess');
%     varsigma{k} = smooth(varsigma{k},0.1,'lowess');
    f = [(mu{k}+2*real(sqrt(varsigma{k})))*scaleVal;flipdim((mu{k}-2*real(sqrt(varsigma{k})))*scaleVal,1)];
    a = fill([Xt; flipdim(Xt,1)], f, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
    a =[ a plot(Xt, mu{k}*scaleVal,'k-')];
    if k>model.nlfPerInt
        c =plot(X{k}+1e-2,y{k}*scaleVal,'k.');
%         d =plot(XGT{k-model.nlfPerInt}, fGT{k-model.nlfPerInt}, 'k--'); 
    end
    minimum = min((mu{k}-2*real(sqrt(varsigma{k})))*scaleVal);
    maximum = max((mu{k}+2*real(sqrt(varsigma{k})))*scaleVal);
    if isfield(model, 'X_u') && ~isempty(model.X_u);
        b = plot(model.X_u, -0.5, 'kx');
        set(b, 'linewidth', 2)
        set(b, 'markersize', 10);
    end
    if k>model.nlfPerInt
        %ylim = [min(minimum,min(y{k}*scaleVal)) max(maximum,max(y{k}*scaleVal))];
        set(c,   'markersize', 0.7*markersize);
%         set(d,   'lineWidth', linewidth);
    else
        %ylim = [min(minimum,min(y{model.nlf+1}*scaleVal)) max(maximum,max(y{model.nlf+1}*scaleVal))];
    end
%        ylabel('y', 'fontsize',fontsize);
%     prop = get(g);
%     poslabel = prop.Position;
%     poslabel(1) = poslabel(1) -0.001*poslabel(1);
%     ylabel('PEV', 'fontsize',fontsize, 'position', poslabel);
%    g = xlabel('Input', 'fontsize',fontsize);
%    prop = get(g);
%    poslabel = prop.Position;
%    poslabel(1) = 0;
%    poslabel(2) = -14;
%    xlabel('Input', 'fontsize',fontsize, 'position', poslabel);
     set(a,   'lineWidth', 2);
     set(gca, 'xlim', xlim, 'ylim', ylim, 'Color', 'none')
    %set(gca, 'position', [0.06 0.08 0.9 0.9])
%    high = get(0, 'screensize');
%    set(gcf, 'position', high)
%    set(gcf, 'PaperPositionMode', 'auto');
    box on
    if saveFigures
        fileName = ['toy3D' upper(model.approx) num2str(k-model.nlf)];
        print('-dpdf', ['../resultsToy1D/' fileName]);
        print('-depsc', ['../resultsToy1D/' fileName], '-loose');
        print('-dpng', ['../resultsToy1D/' fileName]);
    end
end
