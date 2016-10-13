function [params, names] = sdmultiKernExtractParam(kern)

% SDMULTIKERNEXTRACTPARAM Extract parameters from the SDMULTI kernel struc.
% FORMAT
% DESC Extract parameters from the switching dynamical multiple output 
% block kernel structure into a vector of parameters for optimisation. It
% uses the same functions that the multiKernExtractParam.m function but
% also extracts the parameters related to the initial conditions in the
% swicthing dynamical model.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the
% kernel. The vector of transforms is assumed to be empty here, any
% transormation of parameters is assumed to be done in the
% component kernels.
%
% FORMAT
% DESC the same that before but also returns the names of the parameters.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the
% kernel. The vector of transforms is assumed to be empty here, any
% transormation of parameters is assumed to be done in the
% component kernels.
% RETURN names : cell array of strings containing parameter names.
%
% SEEALSO multiKernParamInit, multiKernExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% MODIFICATIONS : Mauricio A. Alvarez, 2010
%
% MODIFICATIONS : Andres F. Lopez Lopera and Mauricio A. Alvarez, 2015

% SDMULTIGP

fhandle = str2func(['sdmulti' kern.comp{1}.type 'KernExtractParam']);
[params, names] =  fhandle(kern);
