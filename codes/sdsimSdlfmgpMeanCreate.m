function meanFunction = sdsimSdlfmgpMeanCreate(q, d, options)

% SDSIMMEANCREATE returns a structure for the switched dynamical SIM mean function.
% FORMAT
% DESC creates the mean function for a multi output GP
% model based in the SDSIM kernel (first order differential equation)
% The outputs of the model are generated according to
%
%       mean_q = (1-cd(t)) B_q/D_q
%
% where mean_q is an output constant corresponding to the mean of the 
% output function q, B_q is basal transcription, D_q is the decay
% constant, and cd(t) = exp(-D_q t)
% RETURN model : the structure for the multigp model
% ARG q : input dimension size.
% ARG d : output dimension size.
% ARG options : contains the options for the MEAN of the multiple output SDSIM model.
%
% SEE ALSO: sdsimKernParamInit, sdsimKernCompute
%
% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio A. Alvarez, 2015

% MULTIGP

if q > 1
  error('SDSIM MEAN FUNCTION only valid for one-D input.')
end

meanFunction.type = 'sdsimSdlfmgp';
meanFunction.nIntervals = length(options.kern.switchingTimes);
meanFunction.nout  = d;
meanFunction.basal = ones(d,1);
meanFunction.decay = ones(d,1);
% meanFunction.switchingTimes = option.kern.switchingTimes;
meanFunction.transforms.index = 1:2*d;
meanFunction.transforms.type = optimiDefaultConstraint('positive');
% Only the parameters of basal rates are counted. The springs are already
% counted in the kernel
meanFunction.nParams = 2*d; 

meanFunction.X = options.X;
meanFunction.switchingTimes = options.kern.switchingTimes;
meanFunction.transforms.index = [meanFunction.transforms.index ...
                                 (2:meanFunction.nIntervals)+meanFunction.transforms.index(end)];
meanFunction.nParams = meanFunction.nParams+meanFunction.nIntervals; 

