function kern = sdmultiKernExpandParam(kern, params)

% SDMULTIKERNEXPANDPARAM Expands parameters into a SDMULTI kernel struct.
% FORMAT
% DESC returns a multiple output block kernel structure filled with the
% parameters in the given vector. This is used as a helper function to
% enable parameters to be optimised in, for example, the NETLAB
% optimisation functions. It has the same functionality than
% multiKernExpandParam.m, but also expands the parameters associated to the
% initial conditions of the SDLFMGP model kernel.
% ARG kern : the kernel structure in which the parameters are to be
% placed.
% ARG param : vector of parameters which are to be placed in the
% kernel structure.
% RETURN kern : kernel structure with the given parameters in the
% relevant locations.
%
% SEEALSO : multiKernParamInit, multiKernExtractParam, kernExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% MODIFICATIONS : Mauricio A. Alvarez, 2010
%
% MODIFICATIONS : Andres F. Lopez Lopera and Mauricio A. Alvarez, 2015

% SDMULTIGP

fhandle = str2func(['sdmulti' kern.comp{1}.type 'KernExpandParam']);
kern =  fhandle(kern, params);

