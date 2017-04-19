function model = sdlfmgpMeanCreate(q , d, options);
% SDLFMGPMEANCREATE creates the mean function for a multi output SDLFM
%
% FORMAT
% DESC returns a structure for the mean function for the multiple output 
% switched dynamical latent force model
% RETURN model : the structure for the mean function of the sdlfmgp model
% ARG q : input dimension size.
% ARG d : output dimension size.
% ARG X : set of training inputs
% ARG y : set of training observations
% ARG options : contains the options for the MEAN of the multi output SDLFM model.
%
% SEE ALSO: sdlfmgpMeanCompute
%
% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio A. Alvarez, 2015

% MULTIGP

fhandle = str2func([options.kernType 'SdlfmgpMeanCreate']);
model   = fhandle(q , d, options);
