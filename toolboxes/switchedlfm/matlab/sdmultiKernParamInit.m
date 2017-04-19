function kern = sdmultiKernParamInit(kern, options)

% SDMULTIKERNPARAMINIT SDMULTI kernel parameter initialisation.
% The swicthing dynamical multiple output block kernel (SDMULTI) is a 
% wrapper kernel designed to represent the situation where there are 
% several Gaussian processes with correlated outputs and initial conditions
% for the system at time zero. We only change some of the features included
% in the multi Kern structure.
% FORMAT
% DESC initialises the switching dynamical multiple output block kernel 
% structure with some default parameters for the initial conditions
% ARG kern : the kernel structure which requires initialisation.
% ARG options : options for the switching dynamical structure.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, multikernParamInit, sdlfmXsdlfmKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2006 
%
% MODIFICATIONS : Mauricio A. Alvarez, 2010
%
% MODIFICATIONS : Andres F. Lopez Lopera and Mauricio A. Alvarez, 2015

% SDMULTIGP


fhandle = str2func(['sdmulti' kern.comp{1}.comp{1}.type 'KernParamInit']);
kern = fhandle(kern, options);
