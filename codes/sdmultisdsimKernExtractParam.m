function [params, names] = sdmultisdsimKernExtractParam(kern)

% SDMULTISDSIMKERNEXTRACTPARAM Extract parameters from the SDMULTI kernel struc.
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

% SDLFMGP

if nargout > 1
  [paramsTemp, namesTemp] = cmpndKernExtractParam(kern);
  namesCovInit = cell(kern.numBlocks);
  
  for i = 1:(1*kern.numPositions)
    for j =1:(1*kern.numPositions)
      namesCovInit{i,j} = ['kInit(' num2str(i) ',' num2str(j) ')'];
    end
  end
  names = {namesTemp{1:kern.nParamsWIC}, namesCovInit{:}};
else
  paramsTemp = cmpndKernExtractParam(kern);
end

% Add the parameters of the covariance for the initial conditions for the
% first interval

params = [paramsTemp(1:kern.nParamsWIC)  kern.LIC(:)'];
