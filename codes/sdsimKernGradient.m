function [g, covGradLocal] = sdsimKernGradient(kern, t, varargin)

% SDSIMKERNGRADIENT Gradient of SDSIM kernel's parameters.
%
% FORMAT
% DESC computes the derivatives as above, but input locations are
% now provided in two matrices associated with rows and columns of
% the kernel matrix. 
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG t1 : the input locations associated with the rows of the
% kernel matrix.
% ARG t2 : the input locations associated with the columns of the
% kernel matrix.
% ARG partial : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The matrix should
% have the same number of rows as t1 and the same number of columns
% as t2 has rows.
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters.
% RETURN covGradLocal : partial derivatives for the initial conditions of
% the first interval.
%
% SEEALSO sdsimKernParamInit, kernGradient, sdsimKernDiagGradient, kernGradX,
% simXsimKernGradient
%
% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio A. Alvarez, 2015

% KERN

if length(varargin)<2
  t2 = t;
else
  t2 = varargin{1};
end

[g1, g2, covGradLocal] = sdsimXsdsimKernGradient(kern, kern, t, t2, varargin{end});

g = real(g1 + g2);
