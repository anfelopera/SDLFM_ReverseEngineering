function K = sdmultiKernCompute(kern, varargin)

% SDMULTIKERNCOMPUTE Compute the SDMULTI kernel given the parameters and X.
% FORMAT
% DESC computes the kernel function for the switching dynamical multiple 
% output block kernel given inputs associated with rows and columns. The
% kernel structure must contain the information about initial conditions.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the switching dynamical multiple 
% output block kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : multiKernParamInit, multiKernCompute.m
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Pei Gao, 2007
%
% MODIFICATIONS : Mauricio A. Alvarez, 2010
%
% MODIFICATIONS : Andres F. Lopez Lopera and Mauricio A. Alvarez, 2015
  
% SDMULTIGP

% Divide the matrix of initial conditions in a matrix of proper dimensions
% for positions and velocities

fhandle = str2func(['sdmulti' kern.comp{1}.type 'KernCompute']);
K =  fhandle(kern, varargin);

