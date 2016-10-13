function K = sdsimXsdsimvKernComputeBlock(simKern1, simKern2, t1, t2, ...
    kyy, i, j, generalConst)

% SDSIMXSDSIMKERNCOMPUTEBLOCK Computes SDSIMxSDSIMV kernel matrix for block i,j
% FORMAT
% DESC computes the kernel matrix for the SDSIMxSDSIMV kernel function in the
% block specified at indeces i,j. It assumes the computation for functions
% that describe positions (system 1) and velocity (system 2)
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
% COPYRIGHT : Andres F. Lopez Lopera and Mauricio A. Alvarez, 2015

% KERN

if nargin<8
    j = i;
    generalConst = [];
end

c1 = sdsimMeanCompute(simKern1(1), t1);
g2 = sdsimvMeanCompute(simKern2(1), t2);

K = kyy*c1*g2.';

if i==j
    for k=1:length(simKern1)
        K  = K + (simvXsimKernCompute(simKern2(k), simKern1(k), t2, t1).');
    end
else    
    if i>j
        PosVel = zeros(1,length(t2));
        for k=1:length(simKern1)
            PosVel = PosVel + simvXsimKernCompute(simKern2(k), simKern1(k), t2, simKern2(k).limit).';
        end
        if isempty(generalConst{i,j})
            K = K + c1*PosVel;
        else
            K = K + (generalConst{i,j}(1,1)*c1)*PosVel;           
        end 
    else
        PosVel = zeros(length(t1),1);
        for k=1:length(simKern1)
            PosVel = PosVel + simXsimKernCompute(simKern1(k), simKern2(k), t1, simKern1(k).limit);            
        end        
         if isempty(generalConst{i,j})
             K = K + PosVel*g2.';
         else
             K = K + PosVel*(generalConst{i,j}(1,1)*g2.');
         end
    end
end
