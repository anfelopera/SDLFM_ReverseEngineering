function K = sdsimXsdsimKernCompute(sdsimKern1, sdsimKern2, t1, t2, covIC, type)

% SDSIMXSDSIMKERNCOMPUTE Compute a cross kernel between two SDSIM kernels.
% FORMAT
% DESC computes cross kernel terms between two SDSIM kernels for
% the multiple output kernel. The SDSIM kernels correspond to Position
% X Position.
% ARG sdsimKern1 : the kernel structure associated with the first SDSIM
% kernel.
% ARG sdsimKern2 : the kernel structure associated with the second SDSIM
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% ARG covIC : covariance for the initial conditions
% ARG type : specifies the type of kerne to be computed
% RETURN K : block of values from kernel matrix.
%
% SEEALSO : sdsimKernParamInit, sdsimKernCompute, sdsimKernParamInit
%
% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio A. Alvarez, 2015

% KERN

if nargin < 6
    type = 'PosPos';
    if nargin < 5
        covIC = t2;
        t2 = t1;
    end
else
    type = 'PosPos';
end

if size(t1, 2) > 1 || size(t2, 2) > 1
    error('Input can only have one column');
end

compareInverseWidth = sdsimKern1.inverseWidth == sdsimKern2.inverseWidth;
if sum(sum(compareInverseWidth))~=(sdsimKern1.nIntervals*sdsimKern1.nlfPerInt)
    error('Kernels cannot be cross combined if they have different inverse widths.')
end
compareSwitchingTimes = sdsimKern1.switchingTimes == sdsimKern2.switchingTimes;
if sum(sum(compareSwitchingTimes))~=sdsimKern1.nIntervals
    error('Kernels cannot be cross combined if they have different switching points.')
end

switch type
    case 'PosPos'
        fhandle = 'sdsimXsdsimKernComputeBlock';
end

fhandle = str2func(fhandle);
    
%Form the basic kernels
simKern1 = struct();
simKern2 = struct();

% Create structures that will make easy the computation of the kernels

spVector = [cumsum(sdsimKern1.switchingTimes) t1(end)+50];

dim1 = zeros(1, sdsimKern1.nIntervals); dim2 = zeros(1, sdsimKern1.nIntervals);

for i=1:sdsimKern1.nIntervals
    for j =1:sdsimKern1.nlfPerInt
        % Create the appropriate set of kernel structures
        simKern1(i,j).decay        = sdsimKern1.decay;     
        simKern1(i,j).inverseWidth = sdsimKern1.inverseWidth(j,i);
        simKern1(i,j).sensitivity  = sdsimKern1.sensitivity(j,i);       
        simKern2(i,j).decay        = sdsimKern2.decay;
        simKern2(i,j).inverseWidth = sdsimKern2.inverseWidth(j,i);
        simKern2(i,j).sensitivity  = sdsimKern2.sensitivity(j,i);
        simKern1(i,j).limit        = spVector(i+1) - spVector(i);
        simKern2(i,j).limit        = spVector(i+1) - spVector(i);
        simKern1(i,j).isNormalised = sdsimKern1.isNormalised;
        simKern2(i,j).isNormalised = sdsimKern2.isNormalised;
        simKern1(i,j).isStationary = sdsimKern1.isStationary;
        simKern2(i,j).isStationary = sdsimKern2.isStationary;
        simKern1(i,j).delay        = sdsimKern1.delay;
        simKern2(i,j).delay        = sdsimKern2.delay;
        simKern1(i,j).isNegativeS  = sdsimKern1.isNegativeS;
        simKern2(i,j).isNegativeS  = sdsimKern2.isNegativeS;        
    end
    newt1 = t1(t1> spVector(i) & t1<spVector(i+1));
    newt2 = t2(t2> spVector(i) & t2<spVector(i+1));
    dim1(i) = length(newt1);
    dim2(i) = length(newt2);
end

if sum(dim1)~=length(t1) || sum(dim2)~=length(t2)
    error('A problem with the dimensions of the switching intervals occured')
end

% Compute some necessary constants

generalConst = sdsimKernComputeConstant(sdsimKern1.nIntervals, ...
    simKern1(1), simKern2(1), spVector);

kyy = zeros(sdsimKern1.nIntervals);
kyy(1,1) = covIC(1,1); 

%%% Compute initial conditions for intervals (1-2), (2-2) and (2-1)

tInit1 = [spVector(1) - spVector(1);spVector(2) - spVector(1)];
tInit2 = [spVector(1) - spVector(1);spVector(2) - spVector(1)];

% Pos -- Pos 
temp = sdsimXsdsimKernComputeBlock(simKern1(1), simKern2(1), ...
    tInit1, tInit2, kyy(1,1), 1, 1, generalConst);
kyy(2,2) = temp(2,2); kyy(1,2) = temp(1,2); kyy(2,1)  = temp(2,1);

tempPosPos = cell(sdsimKern1.nIntervals);

startValOne = 1;
endValOne   = 0;

for i=1:sdsimKern1.nIntervals
    endValOne = endValOne + dim1(i);
    startValThree = 1;
    endValThree = 0;
    for j=1:sdsimKern1.nIntervals
        if i>j
            simKern1Local = simKern1(j,:);
            simKern2Local = simKern2(j,:);
        else
            simKern1Local = simKern1(i,:);
            simKern2Local = simKern2(i,:);
        end            
        endValThree = endValThree + dim2(j);
        % POS -- POS (Kernel and initial positions)
        K(startValOne:endValOne, startValThree:endValThree) = fhandle(simKern1Local, simKern2Local, ...
            t1(startValOne:endValOne) - spVector(i), t2(startValThree:endValThree) - spVector(j), ...
            kyy(i,j), i, j, generalConst);       
        % Time vector initial conditions
        tInit1 = [spVector(i) - spVector(i);spVector(i+1) - spVector(i)];
        tInit2 = [spVector(j) - spVector(j);spVector(j+1) - spVector(j)];        
        tempPosPos{i,j} = sdsimXsdsimKernComputeBlock(simKern1Local, ...
            simKern2Local, tInit1, tInit2, kyy(i,j), i, j, generalConst);
        kyy = organizeIC(kyy, tempPosPos, i, j);
        
        startValThree = endValThree + 1;
    end
    startValOne = endValOne + 1;
end

K = real(K);
