function [grad1, grad2, gsp1, gsp2] = sdsimKernGradientMean(simKern1, ...
    simKern2, t1, t2, kyy, covGrad, typeParam, typeSwitching)

% SDSIMKERNGRADIENTMEAN Gradients of the parameters of mean function cov.
% FORMAT 
% DESC computes the gradients of the parameters that appear in the
% covariance functions formed from the functions that accompany the initial
% conditions.
% ARG simKern1 : structure containing parameters for the system 1
% ARG simKern2 : structure containing parameters for the system 2
% ARG t1 : times at which the system 1 is evaluated
% ARG t2 : times at which the system 2 is evaluated
% ARG kyy : covariance for the initial conditions between position 1 and
% position 2 at block i,j
% ARG covGrad : partial derivatives of the objective function wrt portion
% of the corresponding kernel matrix
% ARG typeParam : specify the mean functions used to compute this part of
% the kernel
% ARG typeParamSwitching : specify the functions used to compute the
% gradients of the swicthing points
% RETURN grad1 : gradients of the parameter of the system 1
% RETURN grad2 : gradients of the parameter of the system 2
% RETURN gsp1 : gradient of the switching point of the system 1
% RETURN gsp2 : gradient of the switching point of the system 2
%
% COPYRIGHT : Andres F. Lopez Lopera and Mauricio A. Alvarez, 2015

% KERN

if nargin<7
    typeParam{1} = 'sdsim';
    typeParam{2} = 'sdsim';
    typeSwitching{1} = 'sdsimv';
    typeSwitching{2} = 'sdsimv';
end

% Computes the derivatives with respect to the parameters in the mean
% function.

fhandle1 = str2func([typeParam{1} 'MeanCompute']);
fhandle2 = str2func([typeParam{2} 'MeanCompute']);

c1 = fhandle1(simKern1, t1);
c2 = fhandle2(simKern2, t2);

fhandle3 = str2func([typeParam{1} 'MeanGradient']);
fhandle4 = str2func([typeParam{2} 'MeanGradient']);

gc1Decay = fhandle3(simKern1, t1);
gc2Decay = fhandle4(simKern2, t2);

% Derivative with respect to Decay1

matGradDecay = kyy*gc1Decay*c2.';

grad1 = (sum(sum(matGradDecay.*covGrad)));

% Derivatives with respect to Decay2

matGradDecay = kyy*c1*gc2Decay.';

grad2 = (sum(sum(matGradDecay.*covGrad)));

% Precomputations for the derivatives of the swicthing points

fhandle5 = str2func([typeSwitching{1} 'MeanCompute']);
fhandle6 = str2func([typeSwitching{2} 'MeanCompute']);

g1 = fhandle5(simKern1, t1);
g2 = fhandle6(simKern2, t2);

% Derivative of switching point 1

matGradSp1 = kyy*g1*c2.';

gsp1 = - sum(sum(matGradSp1.*covGrad));

% Derivative of switching point 2

matGradSp2 = kyy*c1*g2.';

gsp2 = - sum(sum(matGradSp2.*covGrad));

