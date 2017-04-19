function g = sdmultiKernGradient(kern, x, x2, covGrad)

% SDMULTIKERNGRADIENT Gradient of SDMULTI kernel's parameters.
% FORMAT
% DESC computes the gradient of functions with respect to the switching
% dynamical multiple output block
% kernel's parameters. As well as the kernel structure and the
% input positions, the user provides a matrix PARTIAL which gives
% the partial derivatives of the function with respect to the
% relevant elements of the kernel matrix.
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG x : the input locations for which the gradients are being
% computed.
% ARG partial : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The argument takes
% the form of a square matrix of dimension  numData, where numData is
% the number of rows in X.
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters. The ordering of the vector should match
% that provided by the function kernExtractParam.
%
% FORMAT
% DESC computes the derivatives as above, but input locations are
% now provided in two matrices associated with rows and columns of
% the kernel matrix.
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG x1 : the input locations associated with the rows of the
% kernel matrix.
% ARG x2 : the input locations associated with the columns of the
% kernel matrix.
% ARG partial : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The matrix should
% have the same number of rows as X1 and the same number of columns
% as X2 has rows.
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters.
%
% SEEALSO : multiKernParamInit, kernGradient, multiKernDiagGradient, kernGradX
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

fhandle = str2func(['sdmulti' kern.comp{1}.type 'KernGradient']);
if nargin > 3
  g =  fhandle(kern, x, x2, covGrad);
else
  g =  fhandle(kern, x, x2);
end

