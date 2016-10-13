function  model = priorSCreate(q, d, options)
% function  model = meanCreate(q, d, X, y, options)

% PRIORSCREATE creates the prior over the sensitivities for a multi output GP
%
% FORMAT
% DESC returns a structure for the prior over for the multiple output 
% Gaussian process model 
% RETURN model : the structure for the prior over of the multigp model
% ARG q : input dimension size.
% ARG d : output dimension size.
% ARG X : set of training inputs
% ARG y : set of training observations
% ARG options : contains the options for the MEAN of the MULTIGP model.
%
% SEE ALSO: meanCompute
%
% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio A. Alvarez, 2015

% MULTIGP

fhandle = str2func([options.type  'PriorSCreate' ]);
model = fhandle(q , d, options);
model.paramGroups = speye(model.nParams);