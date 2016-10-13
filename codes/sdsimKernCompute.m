function k = sdsimKernCompute(kern, t, t2, covIC)

% SDSIMKERNCOMPUTE Compute the SDSIM kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the swicthing dynamical single input
% motif model kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG t1 : the input matrix associated with the rows of the kernel.
% ARG t2 : the input matrix associated with the columns of the kernel.
% ARG covIC : covariance for the initial conditions
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : sdsimKernParamInit, kernCompute, sdsimKernDiagCompute
%
% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio A. Alvarez, 2015

% KERN


if nargin < 4
  covIC = t2;  
  t2 = t;
end

if size(t, 2) > 1 || size(t2, 2) > 1
  error('Input can only have one column');
end

k = sdsimXsdsimKernCompute(kern, kern, t, t2, covIC);

if nargin < 3;
  k = k + k';
  k = k*0.5;
end

k = real(k); 
