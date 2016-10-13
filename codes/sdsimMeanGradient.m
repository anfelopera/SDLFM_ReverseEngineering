function gradDecay = sdsimMeanGradient(sdsimKern, t)

% SDSIMMEANGRADIENT Gradients wrt parameters of the Position mean SDSIM.
%
% FORMAT
% DESC computes the gradients of the terms $c_d$ that appear in 
% the mean function with respect to decay parameter. 
% ARG sdsimKern : switching dynamical SIM kernel structure with the
% parameters.
% ARG t : input times for which the mean is to be computed.
% RETURN gradDecay : gradient of $c_d$ wrt the Decay
%
% COPYRIGHT : Andres F. Lopez Lopera and Mauricio A. Alvarez, 2015

% SDLFMGP

cd = sdsimMeanCompute(sdsimKern, t);
gradDecay = -t.*cd;  

