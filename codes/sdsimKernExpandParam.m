function kern = sdsimKernExpandParam(kern, params)

% SDLFMKERNEXPANDPARAM Pass parameters from params to SDSIM kernel
% FORMAT
% DESC returns a switching dynamical kernel structure filled with the
% parameters in the given vector. This is used as a helper function to
% enable parameters to be optimised in, for example, the NETLAB
% optimisation functions.
% ARG kern : the kernel structure in which the parameters are to be
% placed.
% ARG param : vector of parameters which are to be placed in the
% kernel structure.
% RETURN kern : kernel structure with the given parameters in the
% relevant locations.
%
% SEEALSO : sdsimKernParamInit, sdsimKernExtractParam, kernExpandParam
%
% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio A. Alvarez, 2015

% KERN

kern.decay = params(1);
kern.inverseWidth = reshape(params(kern.inverseWidthIndx), kern.nlfPerInt, kern.nIntervals);
kern.switchingTimes = params(kern.switchingTimesIndx);
kern.sensitivity = reshape(params(kern.sensitivityIndx), kern.nlfPerInt, kern.nIntervals);
if isfield(kern.options, 'priorS') && kern.options.priorS
    kern.priorSvariance = reshape(params(kern.priorSvarianceIndx), kern.nlfPerInt, kern.nIntervals);
end
