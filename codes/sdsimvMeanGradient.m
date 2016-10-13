function gradDecay = sdsimvMeanGradient(sdsimKern, t)

% SDSIMVMEANGRADIENT Gradients wrt parameters of the velocity mean SDSIM.
% FORMAT
% DESC computes the gradients of the terms $g_d$ that appear in 
% the mean function with respect to the Decay. 
% ARG sdsimKern : switching dynamical SIM kernel structure with the
% parameters.
% ARG t : input times for which the mean is to be computed.
% RETURN gradDecay : gradient of $g_d$ wrt the Decay
%
% COPYRIGHT : Andres F. Lopez Lopera and Mauricio A. Alvarez, 2015

% SDLFMGP

cd = sdsimMeanCompute(sdsimKern, t);
gd = sdsimvMeanCompute(sdsimKern, t);
gradDecay = -t.*gd - cd;
