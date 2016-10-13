function K = sdsimXsdsimKernComputeBlock(simKern1, simKern2, t1, t2, ...
    kyy, i, j, generalConst)

% SDSIMXSDSIMKERNCOMPUTEBLOCK Computes SDSIM kernel matrix for block i,j
% FORMAT
% DESC computes the kernel matrix for the SDSIM kernel function in the
% block specified at indeces i,j. It assumes the computation for functions
% that describe positions (position 1 and position 2).
% ARG simKern1 : structure containing parameters for the system 1
% ARG simKern2 : structure containing parameters for the system 2
% ARG t1 : times at which the system 1 is evaluated
% ARG t2 : times at which the system 2 is evaluated
% ARG kyy : covariance for the initial conditions between position 1 and
% position 2 at block i,j
% ARG i : interval to be evaluated for system 1
% ARG j : interval to be evaluated for system 2
% ARG generalConstant : constants evaluated with sdsimKernComputeConstant.m
% RETURN K : the kernel matrix portion of block i,j
%
% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio A. Alvarez, 2015

% KERN

if nargin<8
    j = i;
    generalConst = [];
end

c1 = sdsimMeanCompute(simKern1(1), t1);
c2 = sdsimMeanCompute(simKern2(1), t2);
K  = kyy*c1*c2.';

if i==j
    for k=1:length(simKern1)
        K  = K + simXsimKernCompute(simKern1(k), simKern2(k), t1, t2);
    end
else
    if i>j
        PosPos = zeros(1,length(t2));
        for k=1:length(simKern1)
            PosPos = PosPos + simXsimKernCompute(simKern1(k), simKern2(k), ...
                simKern2(k).limit, t2);
        end
        if isempty(generalConst{i,j})
            K = K + c1*PosPos;
        else
            K = K + (generalConst{i,j}(1,1)*c1)*PosPos;
        end 
    else
        PosPos = zeros(length(t1),1);
        for k=1:length(simKern1)
            PosPos = PosPos + simXsimKernCompute(simKern1(k), simKern2(k), t1, simKern1(k).limit);
        end        
         if isempty(generalConst{i,j})
             K = K + PosPos*c2.';               
         else
             K = K + PosPos*(generalConst{i,j}(1,1)*c2.');
         end
    end
end

