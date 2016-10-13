function K = sdsimXsdrbfKernComputeBlock(simKern, rbfKern, t1, t2, ...
    i, j, generalConst)

% SDSIMXSDRBFKERNCOMPUTEBLOCK Cross kernel between SDSIM and SDRBF for i,j
% FORMAT
% DESC computes the kernel matrix between a SDSIM kernel function
% and a SDRBF kernel function in the block specified at indeces i,j. It 
% assumes the computation for a function that describe a position.
% ARG simKern : structure containing parameters for the outputs system
% ARG rbfKern : structure containing parameters for the latent system 
% ARG t1 : times at which the system 1 is evaluated
% ARG t2 : times at which the system 2 is evaluated
% ARG i : interval to be evaluated for system 1
% ARG j : interval to be evaluated for system 2
% ARG generalConstant : constants evaluated with sdsimKernComputeConstant.m
% RETURN K : the kernel matrix portion of block i,j
%
% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio A. Alvarez, 2015

% KERN

if nargin<7
    j = i;
    generalConst = [];
end

if i==j
    K  = simXrbfKernCompute(simKern, rbfKern, t1, t2);    
else
    c1 = sdsimMeanCompute(simKern, t1);
    if i>j
        PosRbf = simXrbfKernCompute(simKern, rbfKern, rbfKern.limit, t2);
        if isempty(generalConst{i,j})
            K =  c1*PosRbf;
        else
            K = (generalConst{i,j}(1,1)*c1)*PosRbf;
        end
    end
end
