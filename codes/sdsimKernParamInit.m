function  kern = sdsimKernParamInit(kern)

% SDLFMKERNPARAMINIT SDSIM kernel initialization
% FORMAT
% DESC
% Initializes the switching dynamical single input motif model kernel structure 
% with some initial parameters. The initial parameters are passed through 
% an option in kern.
% RETURN kern : the kernel structure with the default parameters placed in.
% ARG kern : the kernel structure which requires initialisation.
%
% SEEALSO : kernCreate, kernParamInit, lfmKernParamInit
%
% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio A. Alvarez, 2015

% KERN

% Create the basic structure based on the sim kernel.
kern = simKernParamInit(kern);
kern.isNormalised = kern.options.isNormalised;

if isfield(kern.options,'priorS')
    kern.priorS = kern.options.priorS;
else
    kern.priorS = false;
end

% Remove unnecessary fields from the SIM structure
kern = rmfield(kern, {'initVal'});

if isfield(kern, 'options') && isfield(kern.options, 'nIntervals')
    kern.nIntervals = kern.options.nIntervals;
    if isfield(kern.options, 'nlfPerInt')
        kern.nlfPerInt = kern.options.nlfPerInt;
    else
        kern.nlfPerInt = 1;
    end
    kern.inverseWidth = ones(kern.nlfPerInt, kern.nIntervals); % Total number of inverse widths.
    kern.sensitivity = ones(kern.nlfPerInt, kern.nIntervals); % Total number of sensitivities.
    if isfield(kern.options, 'priorS') && kern.options.priorS
        kern.priorSvariance = ones(kern.nlfPerInt, kern.nIntervals);
    end
    % An option for the initialization of the switching times
    if isfield(kern.options, 'switchingTimes')
        if kern.nIntervals == length(kern.options.switchingTimes)
            kern.switchingTimes = kern.options.switchingTimes;
        else
            error('The number of intervals does not match the information of the swicthing time points')
        end
    else
        partition = linspace(-0.1,1, kern.nIntervals + 1);
        kern.switchingTimes = partition(1:end-1);
    end
else
    kern.nIntervals = 1;
    kern.switchingTimes = -0.1;
    warning('SIM:Instead:SDSIM', 'Use the SIM kernel instead.')
end

% Number of parameters computed as:
% decay + inverseWidth (= nlf*intervals)
% + switching times + sensitivities( = nlf*intervals) 
dimParam = [1 kern.nlfPerInt*kern.nIntervals kern.nIntervals ...
    kern.nlfPerInt*kern.nIntervals];
if isfield(kern.options, 'priorS') && kern.options.priorS
    dimParam = [dimParam kern.nlfPerInt*kern.nIntervals];
end

kern.nParams = sum(dimParam);
kern.outputIndx = [1]; % Index of the decay
kern.inverseWidthIndx = dimParam(1)+1:sum(dimParam(1:2));
kern.switchingTimesIndx = sum(dimParam(1:2))+1:sum(dimParam(1:3));
kern.sensitivityIndx = sum(dimParam(1:3))+1:sum(dimParam(1:4));
kern.transforms.index = [1:sum(dimParam(1:2)) (sum(dimParam(1:2))+2):sum(dimParam(1:3))];

if isfield(kern.options, 'priorS') && kern.options.priorS
    kern.priorSvarianceIndx = sum(dimParam(1:4))+1:sum(dimParam);
    kern.transforms.index = [kern.transforms.index sum(dimParam(1:4))+1:sum(dimParam)];
end
kern.transforms.type = optimiDefaultConstraint('positive');
