function meanFunction = sdsimSdlfmgpMeanExpandParam(meanFunction, params)

% SDSIMMEANEXPANDPARAM Extract the parameters of the vector parameter and put 
% them back in the mean function structure for the SDSIM model.
% DESC returns a mean function sim structure filled with the
% parameters in the given vector. This is used as a helper function to
% enable parameters to be optimised in, for example, the NETLAB
% optimisation functions.
% ARG meanFunction : the meanFunction structure in which the parameters are to be
% placed.
% ARG param : vector of parameters which are to be placed in the
% kernel structure.
% RETURN meanFunction : mean function structure with the given parameters in the
% relevant locations.
%
% SEEALSO : sdsimSdlfmgpMeanCreate, sdsimMeanExpandParam, sdsimkernExtractParam
%

% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio A. Alvarez, 2015

% MULTIGP

meanFunction.basal  = params(1:meanFunction.nout)';
meanFunction.decay  = params(meanFunction.nout+1:2*meanFunction.nout)';
meanFunction.switchingTimes = params(2*meanFunction.nout+1:meanFunction.nParams);



