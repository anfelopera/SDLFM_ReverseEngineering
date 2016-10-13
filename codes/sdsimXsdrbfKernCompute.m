function [K, Kc] = sdsimXsdrbfKernCompute(sdsimKern, sdrbfKern, t1, t2, type)

% SDSIMXSDRBFKERNCOMPUTE Cross kernel between a SDSIM and a SDRBF kernels.
% FORMAT
% DESC computes cross kernel terms a SDSIM and a SDRBF kernels for
% the multiple output kernel. The SDSIM kernel corresponds to Position.
% The type of kernel to be computed is specified in 'type'. 
% ARG sdsimKern : the kernel structure associated with the SDSIM
% kernel.
% ARG sdrbfKern : the kernel structure associated with the SDRBF
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% ARG type : specifies the type of kernel to be computed
% RETURN K : block of values from kernel matrix. K is a cell where the
% number of cells is related to the number of latent forces per interval
% RETURN Kc : the kernel matrices grouped in cells. Suitable form for the
% sparse approximations
%
% SEEALSO : sdsimKernParamInit, sdsimKernCompute, sdsimKernParamInit
%
% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio A. Alvarez, 2015

% KERN

if nargin < 5
    type = 'Pos';
    if nargin < 4
        t2 = t1;
    end
else 
    type = 'Pos';
end

if size(t1, 2) > 1 || size(t2, 2) > 1
    error('Input can only have one column');
end

switch type
    case 'Pos'
        fhandle = 'sdsimXsdrbfKernComputeBlock';
end

fhandle = str2func(fhandle);
    
%Form the basic kernels
simKern = struct();
rbfKern = struct();

% Create structures that will make easy the computation of the kernels
%spVector = [sdsimKern.switchingTimes t1(end)+0.1]; % For the last interval include the last point

spVector = [cumsum(sdsimKern.switchingTimes) t1(end)+50];

dim1 = zeros(1, sdsimKern.nIntervals); dim2 = zeros(1, sdrbfKern.nIntervals);


for i =1:sdsimKern.nIntervals
    for j=1:sdsimKern.nlfPerInt
        % Create the appropriate set of kernel structures
        simKern(i,j).decay        = sdsimKern.decay;
        simKern(i,j).inverseWidth = sdsimKern.inverseWidth(j,i);
        simKern(i,j).sensitivity  = sdsimKern.sensitivity(j,i);
        rbfKern(i,j).variance     = sdrbfKern.variance;
        rbfKern(i,j).inverseWidth = sdrbfKern.inverseWidth(j,i);
        simKern(i,j).limit        = spVector(i+1) - spVector(i);
        rbfKern(i,j).limit        = spVector(i+1) - spVector(i);
        simKern(i,j).isNormalised = sdsimKern.isNormalised;
        rbfKern(i,j).isNormalised = sdrbfKern.isNormalised;
        simKern(i,j).isStationary = sdsimKern.isStationary;
        simKern(i,j).isNegativeS   = sdsimKern.isNegativeS;
        simKern(i,j).delay        = sdsimKern.delay;
        simKern(i,j).variance     = sdsimKern.variance;
    end
    newt1   = t1(t1> spVector(i) & t1<spVector(i+1));
    newt2   = t2(t2> spVector(i) & t2<spVector(i+1));
    dim1(i) = length(newt1);
    dim2(i) = length(newt2);
end

if sum(dim1)~=length(t1) || sum(dim2)~=length(t2)
    error('A problem with the dimensions of the switching intervals occured')
end

% Compute some necessary constants

generalConst = sdsimKernComputeConstant(sdsimKern.nIntervals, ...
    simKern(1,1), simKern(1,1), spVector);

K  = cell(sdrbfKern.nlfPerInt,1);
Kc = cell(sdrbfKern.nlfPerInt,1);


for q=1:sdrbfKern.nlfPerInt
    startValOne = 1;
    endValOne   = 0;
    tK = zeros(sum(dim1), sum(dim2));
    for i=1:sdsimKern.nIntervals
        endValOne = endValOne + dim1(i);
        startValThree = 1;
        endValThree = 0;
        for j=1:i
            if i>j
                simKernLocal = simKern(j,q);
                rbfKernLocal = rbfKern(j,q);
            else
                simKernLocal = simKern(i,q);
                rbfKernLocal = rbfKern(i,q);
            end
            endValThree = endValThree + dim2(j);
            tK(startValOne:endValOne, startValThree:endValThree) = fhandle(simKernLocal, ...
                rbfKernLocal, t1(startValOne:endValOne) - spVector(i), ...
                t2(startValThree:endValThree) - spVector(j), i, j, generalConst);
            startValThree = endValThree + 1;
        end
        startValOne = endValOne + 1;
    end
    tK = real(tK);
    K{q} = tK;
end
