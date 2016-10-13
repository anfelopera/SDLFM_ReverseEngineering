function k = sdsimKernDiagComputeBlock(simKern, t, kyy)

% SDSIMKERNCOMPUTEBLOCK Diagonal of a SDSIM kernel matrix for block i
% FORMAT
% DESC computes the Diagonal of a SDSIM kernel matrix for a particular block  
% It assumes the computation for functions that describe positions 
% (position 1 and position 2).
% ARG lfmKern : structure containing parameters for the system
% ARG t : times at which the system is evaluated
% ARG kyy : covariance for the initial conditions between position 1 and
% position 2 at block i,j
% RETURN k : the diagonal of the kernel matrix portion of a block
%
% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio A. Alvarez, 2015

% KERN


c1 = sdsimMeanCompute(simKern(1), t);

k = kyy*(c1.^2);

for i=1:length(simKern)
    k  = k + simKernDiagCompute(simKern(i), t);
end
