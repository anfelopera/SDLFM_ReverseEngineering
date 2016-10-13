function [g1, g2, g3] = sdsimXsdsimKernGradientBlock(simKern1, simKern2, ...
    t1, t2, kyy, i, j, generalConst, generalConstGrad, ...
    covGrad)

% SDSIMXSDSIMKERNGRADIENTBLOCK Gradients of the parameters in block i,j
% FORMAT
% DESC computes the kernel parameters gradients for the SDSIM kernel 
% function in the block specified at indeces i,j. It assumes the 
% computation for functions that describe positions (position 1 and 
% position 2).
% ARG simKern1 : structure containing parameters for the system 1
% ARG simKern2 : structure containing parameters for the system 2
% ARG t1 : times at which the system 1 is evaluated
% ARG t2 : times at which the system 2 is evaluated
% ARG kyy : covariance for the initial conditions between position 1 and
% position 2 at block i,j
% ARG i : interval to be evaluated for system 1
% ARG j : interval to be evaluated for system 2
% ARG generalConstant : constants evaluated with sdsimKernComputeConstant.m
% ARG generalConstGrad : derivatives of the constants computed with
% sdsimKernGradientConstant.m
% ARG covGrad : partial derivatives of the objective function wrt portion
% of the kernel matrix in block i,j
% RETURN g1 : gradients of parameters for the system 1
% RETURN g2 : gradients of parameters for the system 2
% RETURN g3 : gradients of switching points
%
% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio A. Alvarez, 2015.

% KERN

if nargin<8
    j = i;
    generalConst = [];
end

% Compute derivatives of the mean terms with respect to the parameters

[g1Mean, g2Mean, gsp1Mean, gsp2Mean] = sdsimKernGradientMean(simKern1(1), ...
    simKern2(1), t1, t2, kyy, covGrad);

if i==j
    [g1, g2, g3t] = sdsimXsdsimKernGradientBlockIEJ(simKern1, ...
        simKern2, t1, t2, covGrad, g1Mean, g2Mean, gsp1Mean, gsp2Mean);    
    g3(i) = g3t;    

else   
    if i>j
        [g1, g2, g3] = sdsimXsdsimKernGradientBlockIGJ(simKern1, simKern2, ...
            t1, t2, i, j, generalConst, generalConstGrad, covGrad, g1Mean, ...
            g2Mean, gsp1Mean, gsp2Mean); 
   else        
        [g1, g2, g3] = sdsimXsdsimKernGradientBlockILJ(simKern1, simKern2, ...
            t1, t2, i, j, generalConst, generalConstGrad, covGrad, g1Mean, ...
            g2Mean, gsp1Mean, gsp2Mean);
    end
end
