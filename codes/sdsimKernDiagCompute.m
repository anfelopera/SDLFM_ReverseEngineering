function k = sdsimKernDiagCompute(sdsimKern, t, covIC, type)

% SDSIMKERNDIAGCOMPUTE Compute diagonal of a SDSIM kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the switching
% dynamical single input motif model kernel given a design matrix of inputs.
% ARG sdsimKern : the kernel structure for which the matrix is computed.
% ARG t : input data matrix in the form of a design matrix.
% ARG type : specifies the type of diagonal kernel to compute. Options are
% 'PosPos' (default), 'VelVel' and 'AccelAccel'.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : sdlfmKernParamInit, kernDiagCompute, kernCreate
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if nargin < 4
    type = 'PosPos';
end

switch type
    case 'PosPos'
        fhandle = 'sdsimKernDiagComputeBlock';
end

fhandle = str2func(fhandle);


%Form the basic kernels
simKern = struct();

% Create structures that will make easy the computation of the kernels
%spVector = [sdlfmKern.switchingTimes t(end)+0.1]; % For the last interval include the last point

spVector = [cumsum(sdsimKern.switchingTimes) t(end)+50];

dim = zeros(1, sdsimKern.nIntervals);

for i=1:sdsimKern.nIntervals
    for j =1:sdsimKern.nlfPerInt
        % Create the appropriate set of kernel structures
        simKern(i,j).decay = sdsimKern.decay;
        simKern(i,j).inverseWidth = sdsimKern.inverseWidth(i);
        simKern(i,j).sensitivity = sdsimKern.sensitivity(i);
        simKern(i,j).limit = spVector(i+1) - spVector(i);
        simKern(i,j).isNormalised = sdsimKern.isNormalised;
        simKern(i,j).isStationary = sdsimKern.isStationary;
        simKern(i,j).delay = sdsimKern.delay;
        simKern(i,j).isNegativeS = sdsimKern.isNegativeS;        
    end
    newt = t(t> spVector(i) & t<spVector(i+1));
    dim(i) = length(newt);
end


kyy = zeros(sdsimKern.nIntervals,1);
kyy(1,1) = covIC(1,1); 

k = zeros(sum(dim),1);

startValOne = 1;
endValOne   = 0;

for i=1:sdsimKern.nIntervals
    endValOne = endValOne + dim(i);
    k(startValOne:endValOne) = fhandle(simKern(i,:), ...
        t(startValOne:endValOne) - spVector(i), kyy(i));
    tInit = spVector(i+1) - spVector(i);
    kyy(i+1) = sdsimKernDiagComputeBlock(simKern(i,:), tInit, kyy(i));
    startValOne = endValOne + 1;    
end

k = real(k);
