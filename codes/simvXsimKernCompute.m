function [k, sk] = simvXsimKernCompute(simkern1, simkern2, t1, t2)

% SIMVXSIMVKERNCOMPUTE Velocity and regular SIM kernel  
%
% FORMAT
% DESC computes cross kernel terms between SIMV and SIM kernels for
% the multiple output kernel.
% ARG simKern1 : the kernel structure associated with the first SIMV
% kernel.
% ARG simKern2 : the kernel structure associated with the second SIM
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% RETURN K : block of values from kernel matrix.
%
% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio A. Alvarez, 2015

% KERN

if nargin < 4
    t2 = t1;
end
if size(t1, 2) > 1 || size(t2, 2) > 1
    error('Input can only have one column');
end

if simkern1.inverseWidth ~= simkern2.inverseWidth
    error('Kernels cannot be cross combined if they have different inverse widths.')
end

sigma = sqrt(2/simkern1.inverseWidth);

if (simkern1.isStationary == false)
    [~,~,~,~,dh1_dt] = simComputeH(t1, t2, simkern1.decay, simkern2.decay, simkern1.delay, simkern2.delay, sigma);
else
    [~,~,~,~,dh1_dt] = simComputeHStat(t1, t2, simkern1.decay, simkern2.decay, simkern1.delay, simkern2.delay, sigma);
end
if nargin < 4
  sk = 0.5 * (dh1_dt + dh1_dt');
else
  if (simkern2.isStationary == false)
      [~,~,~,~,~,dh2_dt] = simComputeH(t2, t1, simkern2.decay, simkern1.decay, simkern2.delay, simkern1.delay, sigma);
  else
      [~,~,~,~,~,dh2_dt] = simComputeHStat(t2, t1, simkern2.decay, simkern1.decay, simkern2.delay, simkern1.delay, sigma);
  end
  sk = 0.5 * (dh1_dt + dh2_dt');
end
if isfield(simkern1, 'isNormalised') && (simkern1.isNormalised == true)
    k = sk;
else
    k = sqrt(pi)*sigma*sk;
end

if isfield(simkern1, 'isNegativeS') && (simkern1.isNegativeS == true)
    k = (simkern1.sensitivity*simkern2.sensitivity)*k;
else
    k = simkern1.variance*k;
end

if isfield(simkern1, 'gaussianInitial') && simkern1.gaussianInitial,
  dim1 = size(t1, 1);
  dim2 = size(t2, 1);
  t1Mat = t1(:, ones(1, dim2));
  t2Mat = t2(:, ones(1, dim1))';

  k = k + simkern1.initialVariance * exp(- kern.decay * (t1Mat + t2Mat));
end
