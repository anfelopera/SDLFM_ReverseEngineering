function covGradLocal = sdsimKernMeanCovPartial(simKern1, simKern2, t1, ...
    t2, covGrad, typeParam)

% SDSIMKERNMEANCOVPARTIAL Helper function for derivatives in SDSIM kernel
% FORMAT
% DESC computes partial derivatives for the total gradinets of the
% parameters in the SDSIM kernel
% ARG simKern1 : structure containing parameters system 1
% ARG simKern2 : structure containing parameters system 2
% ARG t1 : time points for the evaluation of system 1
% ARG t2 : time points for the evaluation of system 2
% ARG covGrad : portion of the partial derivative with respect the
% objective function that it's being optimised.
% ARG typeParam : associated to the type of kernel being computed.
% RETURN covGradLocal : the partial derivatives of the parameters

% KERN

if nargin < 6
    typeParam = {'sdsim', 'sdsim'};
end

fhandle1 = str2func([typeParam{1} 'MeanCompute']);
fhandle2 = str2func([typeParam{2} 'MeanCompute']);

c1 = fhandle1(simKern1, t1);
c2 = fhandle2(simKern2, t2);

covGradLocal(1,1) = sum(sum((c1*c2').*covGrad));
