function [g1, g2, g3] = sdsimXsdsimKernGradientBlockILJ(simKern1, simKern2, ...
    t1, t2, i, j, generalConst, generalConstGrad, covGrad, g1Mean, ...
    g2Mean, gsp1Mean, gsp2Mean, typeMeanParam, typeMeanSwitching, ...
    typeKernParam1, typeKernSwitching1)

% SDSIMXSDSIMKERNGRADIENTBLOCKILJ 
% FORMAT
% DESC computes the gradients of the parameters for system 1 and system 2
% when i is lower than j.
% ARG simKern1 : structure containing parameters for the system 1
% ARG simKern2 : structure containing parameters for the system 2
% ARG t1 : times at which the system 1 is evaluated
% ARG t2 : times at which the system 2 is evaluated
% ARG i : interval to be evaluated for system 1
% ARG j : interval to be evaluated for system 2
% ARG generalConstant : constants evaluated with sdsimKernComputeConstant.m
% ARG generalConstGrad : derivatives of the constants computed with
% sdsimKernGradientConstant.m
% ARG covGrad : partial derivatives of the objective function wrt portion
% of the corresponding kernel matrix
% ARG g1Mean : gradients of the parameter of the system 1 obtained from the
% part of the function that uses the funcitons accommpanying the initial
% conditions.
% ARG g2Mean : gradients of the parameter of the system 2 obtained from the
% part of the function that uses the funcitons accommpanying the initial
% conditions.
% ARG gsp1Mean : gradient of the switching point of the system 1 obtained 
% from the part of the function that uses the funcitons accommpanying the 
% initial conditions.
% ARG gsp2Mean : gradient of the switching point of the system 2 obtained 
% from the part of the function that uses the funcitons accommpanying the 
% initial conditions.
% ARG typeMeanParam : specify the mean functions used to compute this part 
% of the kernel
% ARG typeMeanSwitching : specify the functions used to compute the
% gradients of the swicthing points in part thst uses mean functions
% ARG typeKernParam1 : specify the first kernel function used to compute 
% this part of the kernel
% ARG typeKernSwitching1 : specify the functions used to compute the
% gradients of the swicthing points in both sides of the kernel function 1
% RETURN g1 : gradients of parameters for the system 1
% RETURN g2 : gradients of parameters for the system 2
% RETURN g3 : gradients of switching points
%
% COPYRIGHT : Mauricio A. Alvarez, 2010.

% KERN

if nargin < 14
    typeMeanParam = 'sdsim';
    typeMeanSwitching = 'sdsimv';
    typeKernParam1 = 'simXsim';
    typeKernSwitching1= {'simvXsim', 'simvXsim'};
end

g3 = [];

fhandleMeanParam = str2func([typeMeanParam 'MeanCompute']);
fhandleMeanGradParam = str2func([typeMeanParam 'MeanGradient']);
fhandleMeanGradSwitching = str2func([typeMeanSwitching 'MeanCompute']);
fhandleKernPosPos = str2func([typeKernParam1 'KernCompute']);
fhandleKernGradPosPos = str2func([typeKernParam1 'KernGradient']);
fhandleKernGradSwitchingPosPos1 = str2func([typeKernSwitching1{1} 'KernCompute']);
fhandleKernGradSwitchingPosPos2 = str2func([typeKernSwitching1{2} 'KernCompute']);


c2 = fhandleMeanParam(simKern2(1), t2);

if isempty(generalConst{i,j})
    coeffPosPos = c2;
else
    coeffPosPos = generalConst{i,j}(1,1)*c2;
end
g1Kern = zeros(length(simKern1), 3);
g2Kern = zeros(length(simKern1), 3);
PosPos = zeros(length(t1),1);
for k=1:length(simKern1)
    PosPos = PosPos + fhandleKernPosPos(simKern1(k), simKern2(k), t1,...
        simKern1(k).limit);
    [g1KLocal, g2KLocal] = fhandleKernGradPosPos(simKern1(k), simKern2(k), ...
        t1, simKern1(k).limit, covGrad, coeffPosPos.');
    g1Kern(k,:) = g1Kern(k, :) + g1KLocal;
    g2Kern(k,:) = g2Kern(k, :) + g2KLocal;
end
gcDecay = fhandleMeanGradParam(simKern2(1), t2);
gradDecay = 1;
g2Pos = fhandleMeanGradSwitching(simKern2(1), t2);
gsp1Kern = zeros(1, length(simKern1));
gsp2Kern = zeros(1, length(simKern1));
% This means that are only two switching points involved, t1
% and t0 appear in k, like k_{ff}(t1-t0,t'-t0) and
% k_{mf}(t1-t0,t'-t0). The other points appear in c1(t - t1)
for k=1:length(simKern1)
    temp = fhandleKernGradSwitchingPosPos1(simKern1(k), simKern2(k), t1, simKern1(k).limit);
    gsp1Kern(k) = gsp1Kern(k) - sum(sum((temp*coeffPosPos.').*covGrad));
    temp = fhandleKernGradSwitchingPosPos2(simKern2(k), simKern1(k), simKern1(k).limit, t1).';
    gsp1Kern(k) = gsp1Kern(k) - sum(sum((temp*coeffPosPos.').*covGrad));
    gsp2Kern(k) = gsp2Kern(k) + sum(sum((temp*coeffPosPos.').*covGrad));
end
if isempty(generalConst{i,j})
    matGradDecay = PosPos*gcDecay.';
    gCoeff = gradDecay*sum(sum(matGradDecay.*covGrad));
    matGradSp2 = PosPos*g2Pos.';
    gsp2 = - sum(sum(matGradSp2.*covGrad));
    g3(i) = gsp1Mean + sum(gsp1Kern);                % switching points 1
    g3(j) = gsp2Mean + sum(gsp2Kern) + gsp2;         % switching points 2
else
    constGradDecay = generalConstGrad{1}{i,j};
    constVal = generalConst{i,j};
    % Derivative wrt parameters
    matGradDecay = PosPos*((constVal(1,1)*gcDecay.' + constGradDecay(1,1)*c2.'));
    gCoeff = gradDecay*sum(sum(matGradDecay.*covGrad));
    % Derivative wrt switching points
    % Firts, notice that gsp2Kern doesn't correspond to the switching point
    % of the interval j (or say 2). It's just the derivative of the inner
    % switching point that lead to the actual switching point for the
    % interval j. So we rename it first.
    gspInt = sum(gsp2Kern);
    % Compute the derivative of the swicthing point j, not appearing in the
    % constant
    matGradSp2 = PosPos*(constVal(1,1)*g2Pos.');
    gsp2 = - sum(sum(matGradSp2.*covGrad));
    % Compute the derivatives for all the switching points apearing in the
    % constant
    constGradSPoint = generalConstGrad{2}{i,j};
    numberSP = size(constGradSPoint,2);
    gspInBetween = zeros(1, numberSP);
    for k=1:numberSP
        temp = PosPos*(constGradSPoint(1,k)*c2.');
        gspInBetween(k) = sum(sum(temp.*covGrad));
    end
    % Assign derivatives wrt all other switching points, with a correction
    % for the derivative of the innermost switching point in the constant
    gspInBetween(end) =  gspInBetween(end) + gspInt;
    g3(i) = gsp1Mean + sum(gsp1Kern);           % switching point 1
    g3(j) = gsp2Mean + gsp2 + gspInBetween(1);  % switching point 2
    g3(i+1:j-1) = fliplr(gspInBetween(2:end));
end
% Assign derivatives wrt first system
g1{1} = g1Mean + sum(g1Kern(:,1), 1);             % decay 1
g1{2} = g1Kern(:,2).';                            % inverse widths
g1{3} = g1Kern(:,3).';                            % sensitivities 1
% Assign derivatives wrt second system
g2{1} = g2Mean + sum(g2Kern(:,1), 1) + gCoeff;    % decay 2
g2{2} = g2Kern(:,2).';                            % inverse widths
g2{3} = g2Kern(:,3).';                            % sensitivities 2
