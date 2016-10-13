function model = sdsimSdlfmgpFixParam(model)

% SDLFMSDLFMGPFIXPARAM Fix parameters for a sdlfmgp gp with sdsim kernel
% FORMAT
% DESC Fix the parameters for a sdlfmgp model that uses SDSIM kernel.
% RETURN model : model with fixed parameters included
% ARG model    : model before fixing the parameters
%
% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio A. Alvarez, 2015

% SDLFMGP

count = 1;
model.fix(count).index = [];
model.fix(count).value = [];
