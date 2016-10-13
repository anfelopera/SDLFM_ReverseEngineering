function k = sdmultiKernDiagCompute(kern, x)

% SDMULTIKERNDIAGCOMPUTE Compute diagonal of SDMULTI kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the switching
% dynamical multiple output block kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : sdmultiKernParamInit, kernDiagCompute, sdmultiKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Pei Gao, 2007

% MODIFICATIONS : Mauricio Alvarez, 2008, 2010

% MODIFICATIONS : Andres F. Lopez Lopera and Mauricio A. Alvarez, 2015

% SDMUTLIGP

% It actually calls the same function that already exists for multiKern

% Divide the matrix of initial conditions in a matrix of proper dimensions
% for positions and velocities

fhandle = str2func(['sdmulti' kern.comp{1}.type 'KernDiagCompute']);
k =  fhandle(kern, x);
