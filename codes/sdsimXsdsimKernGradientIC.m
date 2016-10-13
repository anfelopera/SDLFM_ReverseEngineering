function [g1, g2, g3, g4] = sdsimXsdsimKernGradientIC(simKern1, simKern2, ...
    t1, t2, gkyy1, gkyy2, gkyy3, gkyy4, covGrad, typeParam)

% SDSIMXSDSIMKERNGRADIENTIC Computes partial derivative for init. const.
%
% COPYRIGHT : Andres f. Lopez Lopera and Mauricio A. Alvarez, 2015

% KERN

if nargin < 10
    typeParam = {'sdsim', 'sdsim'};    
end

g1 = cell(1,3);
g2 = cell(1,3);

covGradLocal = sdsimKernMeanCovPartial(simKern1(1), simKern2(1), t1, t2, covGrad, ...
    typeParam);

% Parameters

for i=1:3,
    g1{i} = gkyy1{i}*covGradLocal(1,1);
    g2{i} = gkyy2{i}*covGradLocal(1,1);
end

% Switching points

g3 = gkyy3*covGradLocal(1,1);

% First covariance

g4 = gkyy4*covGradLocal(1,1);