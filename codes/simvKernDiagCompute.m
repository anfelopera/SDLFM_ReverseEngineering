function k = simvKernDiagCompute(kern, t)

% SIMVKERNDIAGCOMPUTE Compute diagonal of SIMV kernel.
%
%	Description:
%
%	K = SIMVKERNDIAGCOMPUTE(KERN, T) computes the diagonal of the kernel
%	matrix for the derivative of the single input motif kernel given a design
%   matrix of
%	inputs.
%	 Returns:
%	  K - a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the matrix is computed.
%	  T - input data matrix in the form of a design matrix.
%	
%
%	See also
%	KERNDIAGCOMPUTE, KERNCREATE, SIMKERNCOMPUTE

% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio A. Alvarez, 2015


if size(t, 2) > 1 
  error('Input can only have one column');
end

sigma = sqrt(2/kern.inverseWidth);
t = t - kern.delay;
halfSigmaD = 0.5*sigma*kern.decay;

if (kern.isStationary == false)
    lnPart1 = lnDiffErfs(halfSigmaD, halfSigmaD - t/sigma);
    dh_dt = +exp(halfSigmaD*halfSigmaD-(t/sigma+halfSigmaD).*2)...
        +exp(halfSigmaD*halfSigmaD-2*kern.decay*t).*...
        (exp(lnPart1+log(2*kern.decay)) - exp(-(t/sigma-halfSigmaD).*2) );
else
    dh_dt = zeros(size(t));
end

if isfield(kern, 'isNegativeS') && (kern.isNegativeS == true)
    k = (kern.sensitivity*kern.sensitivity)*dh_dt/(2*kern.decay);
else
    k = kern.variance*dh_dt/(2*kern.decay);
end

if ~isfield(kern, 'isNormalised') || (kern.isNormalised == false)
    k = sqrt(pi)*sigma*k;
end

if isfield(kern, 'gaussianInitial') && kern.gaussianInitial,
  k = k + kern.initialVariance*exp(-2*kern.decay*t);
end
