function meanFunction = lfmMeanExpandParam(meanFunction, params)

% LFMMEANEXPANDPARAM Extract the parameters of the vector parameter and put 
% them back in the mean function structure for the LFM model.
% DESC returns a mean function lfm structure filled with the
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
% SEEALSO : lfmMeanCreate, lfmMeanExtractParam, kernExpandParam
%
% COPYRIGHT : Mauricio Alvarez and Neil D. Lawrence, 2008

% MULTIGP

meanFunction.basal  = params(1:meanFunction.nParams/2)';
meanFunction.spring = params(meanFunction.nParams/2+1:meanFunction.nParams)';

