function generalConst = sdsimKernComputeConstant(nIntervals, ...
    simKern1, simKern2, spVector)

% SDSIMKERNCOMPUTECONSTANT Compute constants for the SDSIM kernel
% FORMAT
% DESC computes necessary constants in order to compute the SDSIM kernel
% matrix.
% ARG nIntervals : number of switching intervals in the kernel
% ARG simKern1 : structure containing the parameters of system 1
% ARG simkern2 : structure containing the parameters of system 2
% ARG spVector : vector containing the switching time values
% RETURN generalConstant : a cell containing the necessary constants for
% computing the kernel in the switching intervals.
%
% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio A. Alvarez, 2015

% KERN

generalConst = cell(nIntervals);

for i=1:nIntervals
    for j =1:i-1
        if i - j>=2
            c1 = sdsimMeanCompute(simKern1, spVector(i) - spVector(i-1));
            c2 = sdsimMeanCompute(simKern2, spVector(i) - spVector(i-1));
            g1 = sdsimvMeanCompute(simKern1, spVector(i) - spVector(i-1));
            g2 = sdsimvMeanCompute(simKern2, spVector(i) - spVector(i-1));
            if i - j == 2
                generalConst{i,j}(1,1) = c1;
                generalConst{i,j}(2,1) = g1;
                generalConst{j,i}(1,1) = c2;
                generalConst{j,i}(2,1) = g2;
            else
                generalConst{i,j}(1,1) = c1*generalConst{i-1,j}(1,1);
                generalConst{j,i}(1,1) = c2*generalConst{j,i-1}(1,1);
                generalConst{i,j}(2,1) = g1*generalConst{i-1,j}(1,1);
                generalConst{j,i}(2,1) = g2*generalConst{j,i-1}(1,1);
            end
        end
    end
end
