function [g1, g2, covGradLocal] = sdsimXsdsimKernGradient(sdsimKern1, ...
    sdsimKern2, t1, t2, covGrad, covIC, type)

% SDSIMXSDSIMKERNGRADIENT Gradients of cross kernel between 2 SDSIM kernels.
% FORMAT
% DESC computes a cross gradient for a cross kernel between two switching 
% dynamical SIM kernels for the multiple output kernel. SDSIM kernels can 
% correspond to Position X Position (default). 
% The type of kernel for which the gradients re obtained is specified in 'type'. 
% ARG sdsimKern1 : the kernel structure associated with the first SDSIM
% kernel.
% ARG sdsimKern2 : the kernel structure associated with the second SDSIM
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% ARG covIC : covariance for the initial conditions
% ARG type : specifies the type of kernel to be computed
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see simKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see simKernExtractParam.
% RETURN covGrad : partial covariances for kyy(t_0, t_0)
%
% SEEALSO : sdsimKernParamInit, sdsimKernCompute, sdsimKernParamInit
%
% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio A. Alvarez, 2015

% KERN

if nargin == 5
    covIC = covGrad;
    covGrad = t2;
    t2 = t1;
    type = 'PosPos';
elseif nargin == 6
    type = 'PosPos';
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
        fhandle1 = 'sdsimXsdsimKernGradientBlock';
        typeIC   = {'sdsim', 'sdsim'};
end

fhandle = str2func(fhandle1);

%Form the basic kernels
simKern1 = struct();
simKern2 = struct();

% Create structures that will make easy the computation of the kernels

spVector = [cumsum(sdsimKern1.switchingTimes) t1(end)+100];

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
        if isfield(sdsimKern1,'priorS') && isfield(sdsimKern2,'priorS')...
                && sdsimKern1.priorS && sdsimKern2.priorS
            simKern1(i,j).priorS          = sdsimKern1.priorS;        
            simKern2(i,j).priorS          = sdsimKern2.priorS;          
            simKern1(i,j).priorSvariance  = sdsimKern1.priorSvariance(j,i);                
            simKern2(i,j).priorSvariance  = sdsimKern2.priorSvariance(j,i);        
        end
    end
    newt1 = t1(t1> spVector(i) & t1<spVector(i+1));
    newt2 = t2(t2> spVector(i) & t2<spVector(i+1));
    dim1(i) = length(newt1);
    dim2(i) = length(newt2);
end

% Compute some necessary constants and their gradients

[generalConstGrad, generalConst] = sdsimKernGradientConstant(sdsimKern1.nIntervals, ...
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
g1PP = cell(sdsimKern1.nIntervals); g2PP = cell(sdsimKern1.nIntervals);
g3PP = cell(sdsimKern1.nIntervals); g4PP = cell(sdsimKern1.nIntervals);
gkyy1 = cell(sdsimKern1.nIntervals); gkyy2 = cell(sdsimKern1.nIntervals);
gkyy3 = cell(sdsimKern1.nIntervals); gkyy4 = cell(sdsimKern1.nIntervals);

startValOne = 1;
endValOne   = 0;


% Initialization of the vector of gradients
g1 = zeros(1, sdsimKern1.nParams);
g2 = zeros(1, sdsimKern1.nParams);
gradS1 = zeros(sdsimKern1.nlfPerInt, sdsimKern1.nIntervals);
gradS2 = zeros(sdsimKern1.nlfPerInt, sdsimKern1.nIntervals);
gradIW1 = zeros(sdsimKern1.nlfPerInt, sdsimKern1.nIntervals);
gradIW2 = zeros(sdsimKern1.nlfPerInt, sdsimKern1.nIntervals);
gradSP = zeros(1, sdsimKern1.nIntervals);


%%% Precomputations

% Needed to compute derivatives with respect to the covariances of the
% initial conditions

covGradLocal = sdsimKernMeanCovPartial(simKern1(1), simKern2(1), ...
    t1(1:dim1(1)) - spVector(1), t2(1:dim2(1)) - spVector(1), ...
    covGrad(1:dim1(1), 1:dim2(1)), typeIC);

% Make initialitions for some derivatives
for k=1:3
    gkyy1{1,1}{k} = 0; gkyy2{1,1}{k} = 0;
end

gkyy3{1,1} = 0;
gkyy4{1,1} = zeros(2);


for i=1:sdsimKern1.nIntervals
    endValOne = endValOne + dim1(i);
    startValThree = 1;
    endValThree = 0;
    for j=1:sdsimKern1.nIntervals
        if i>j
            simKern1Local = simKern1(j,:);
            simKern2Local = simKern2(j,:);
            lowest = j;
            biggest = i;
        else
            simKern1Local = simKern1(i,:);
            simKern2Local = simKern2(i,:);
            lowest = i;
            biggest = j;
        end
        endValThree = endValThree + dim2(j);
        % POS -- POS (Kernel and initial positions)
        [g1Local, g2Local, g3Local] = fhandle(simKern1Local, simKern2Local, ...
             t1(startValOne:endValOne) - spVector(i), t2(startValThree:endValThree) - spVector(j), ...
            kyy(i,j), i, j, generalConst, generalConstGrad, ...
            covGrad(startValOne:endValOne, startValThree:endValThree));
        % Time vector initial conditions
        tInit1 = [spVector(i) - spVector(i);spVector(i+1) - spVector(i)];
        tInit2 = [spVector(j) - spVector(j);spVector(j+1) - spVector(j)];        
        tempPosPos{i,j} = sdsimXsdsimKernComputeBlock(simKern1Local, ...
            simKern2Local, tInit1, tInit2, kyy(i,j), i, j, generalConst);
        % Gradients         
        kyy = organizeIC(kyy, tempPosPos, i, j);
        [g1PP{i,j}, g2PP{i,j}, g3PP{i,j}, g4PP{i,j}] = sdsimXsdsimKernGradientICBlock(simKern1Local, ...
            simKern2Local, tInit1, tInit2, kyy(i,j), i, j, ...
            generalConst, generalConstGrad,gkyy1, gkyy2, gkyy3, gkyy4);
        % Organize derivatives
        [gkyy1, gkyy2, gkyy3, gkyy4] = organizeDerivatives(gkyy1, gkyy2, gkyy3, gkyy4, ...
            g1PP, g2PP, g3PP , g4PP, i,j);        
        % Be aware that these derivatives are from the parameters in
        % the interval before. This is important for inversewidths,
        % switching points and sensitivities. For the other parameters
        % it is not so important.        
        [g1LIC, g2LIC, g3LIC, g4LIC] = sdsimXsdsimKernGradientIC(simKern1Local(1), simKern2Local(1), ...
            t1(startValOne:endValOne) - spVector(i), t2(startValThree:endValThree) - spVector(j), ...
            gkyy1{i,j}, gkyy2{i,j}, gkyy3{i,j}, gkyy4{i,j}, covGrad(startValOne:endValOne, ...
            startValThree:endValThree), typeIC);
        % Organize the normal derivatives
        g1(sdsimKern1.outputIndx) = g1(sdsimKern1.outputIndx) + g1Local{1} + g1LIC{1};
        g2(sdsimKern1.outputIndx) = g2(sdsimKern1.outputIndx) + g2Local{1} + g2LIC{1};
        if i>=2 && j>=2
            localIW = reshape(g1LIC{2},sdsimKern1.nlfPerInt, lowest-1);
            gradIW1(:,1:lowest) = gradIW1(:,1:lowest) + [localIW  g1Local{2}'];
            localIW = reshape(g2LIC{2},sdsimKern1.nlfPerInt, lowest-1);
            gradIW2(:,1:lowest) = gradIW2(:,1:lowest) + [localIW g2Local{2}'];
            localS = reshape(g1LIC{3},sdsimKern1.nlfPerInt, lowest-1);
            gradS1(:,1:lowest) = gradS1(:,1:lowest) + [localS g1Local{3}'];
            localS = reshape(g2LIC{3},sdsimKern1.nlfPerInt, lowest-1);
            gradS2(:,1:lowest) = gradS2(:,1:lowest) + [localS g2Local{3}'];
        else
            gradIW1(:,lowest) = gradIW1(:,lowest) + g1Local{2}';
            gradIW2(:,lowest) = gradIW2(:,lowest) + g2Local{2}';
            gradS1(:,lowest) = gradS1(:,lowest) + g1Local{3}';
            gradS2(:,lowest) = gradS2(:,lowest) + g2Local{3}';
        end
        gradSP(1:biggest) = gradSP(1:biggest) + g3Local + g3LIC;
        covGradLocal = covGradLocal + g4LIC;
        startValThree = endValThree + 1;
    end
    startValOne = endValOne + 1;
end

g1(sdsimKern1.inverseWidthIndx) = gradIW1(:)' + gradIW2(:)';
g2(sdsimKern1.inverseWidthIndx) = zeros(1, sdsimKern1.nlfPerInt*sdsimKern1.nIntervals);
tempGradSP = fliplr(gradSP);
tempGradSP = cumsum(tempGradSP);
g1(sdsimKern1.switchingTimesIndx) = fliplr(tempGradSP); 
g2(sdsimKern1.switchingTimesIndx) = zeros(1, sdsimKern1.nIntervals);
g1(sdsimKern1.sensitivityIndx) = gradS1(:)';
g2(sdsimKern1.sensitivityIndx) = gradS2(:)';

if isfield(sdsimKern1, 'priorS') && sdsimKern1.priorS
    nout = length(g1LIC);
    z = 1/(nout*sdsimKern1.options.nlfPerInt); 
    
    g1(sdsimKern1.sensitivityIndx) = g1(sdsimKern1.sensitivityIndx)'...
         - z*sdsimKern1.sensitivity(:)./sdsimKern1.priorSvariance(:);
    g2(sdsimKern1.sensitivityIndx) = g2(sdsimKern1.sensitivityIndx)'...
         - z*sdsimKern2.sensitivity(:)./sdsimKern2.priorSvariance(:);
    
    g1(sdsimKern1.priorSvarianceIndx) = z*(- 0.5*(1./sdsimKern1.priorSvariance(:))...
        + 0.5*(sdsimKern1.sensitivity(:)./sdsimKern1.priorSvariance(:)).^2);
    g2(sdsimKern1.priorSvarianceIndx) = z*(- 0.5*(1./sdsimKern2.priorSvariance(:))...
        + 0.5*(sdsimKern2.sensitivity(:)./sdsimKern2.priorSvariance(:)).^2);
end

function [gkyy1, gkyy2, gkyy3, gkyy4] = organizeDerivatives(gkyy1, gkyy2, gkyy3, gkyy4, ...
    g1PP, g2PP, g3PP , g4PP, i,j)

if i==1 && j==1
    gkyy1{i+1, j+1} = g1PP{i,j}{2,2}; gkyy2{i+1, j+1} = g2PP{i,j}{2,2};
    gkyy3{i+1, j+1} = g3PP{i,j}{2,2}; gkyy4{i+1, j+1} = g4PP{i,j}{2,2};
    gkyy1{i,j+1} = g1PP{i,j}{1,2}; gkyy2{i,j+1} = g2PP{i,j}{1,2};
    gkyy3{i,j+1} = g3PP{i,j}{1,2}; gkyy4{i,j+1} = g4PP{i,j}{1,2};
    gkyy1{i+1,j} = g1PP{i,j}{2,1}; gkyy2{i+1,j} = g2PP{i,j}{2,1};
    gkyy3{i+1,j} = g3PP{i,j}{2,1}; gkyy4{i+1,j} = g4PP{i,j}{2,1}; 
end
if i==1 && j~=1
    gkyy1{i+1, j+1} = g1PP{i,j}{2,2}; gkyy2{i+1, j+1} = g2PP{i,j}{2,2};
    gkyy3{i+1, j+1} = g3PP{i,j}{2,2}; gkyy4{i+1, j+1} = g4PP{i,j}{2,2};
    gkyy1{i,j+1} = g1PP{i,j}{1,2}; gkyy2{i,j+1} = g2PP{i,j}{1,2};
    gkyy3{i,j+1} = g3PP{i,j}{1,2}; gkyy4{i,j+1} = g4PP{i,j}{1,2};
end
if i~=1 && j==1
    gkyy1{i+1, j+1} = g1PP{i,j}{2,2}; gkyy2{i+1, j+1} = g2PP{i,j}{2,2};
    gkyy3{i+1, j+1} = g3PP{i,j}{2,2}; gkyy4{i+1, j+1} = g4PP{i,j}{2,2};
    gkyy1{i+1,j} = g1PP{i,j}{2,1}; gkyy2{i+1,j} = g2PP{i,j}{2,1};
    gkyy3{i+1,j} = g3PP{i,j}{2,1}; gkyy4{i+1,j} = g4PP{i,j}{2,1}; 
end
if i~=1 && j~=1
    gkyy1{i+1, j+1} = g1PP{i,j}{2,2}; gkyy2{i+1, j+1} = g2PP{i,j}{2,2};
    gkyy3{i+1, j+1} = g3PP{i,j}{2,2}; gkyy4{i+1, j+1} = g4PP{i,j}{2,2};
end
