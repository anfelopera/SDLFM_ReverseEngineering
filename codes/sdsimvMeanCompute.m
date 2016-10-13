function term  = sdsimvMeanCompute(sdsimKern, t)

% SDSIMVMEANCOMPUTE Velocity mean for the switching dynamical SIM model.
% Computes the terms $g_d$ that appear in the mean function associated 
% with the switching dynamical SIM model. If the mean function 
% is mu(t), then
%
%     mu(t) = g_d(t)y_d(t_0),
%
% where $y_d(t_0)$ is the initial condition associated to the position.
% 
% FORMAT
% DESC
% ARG sdsimKern : switching dynamical SIM kernel structure with the
% parameters.
% ARG t : input times for which the mean is to be computed.
% RETURN term : the value of $g_d$.
%
% FORMAT
% DESC
% Computes the terms that appear in the mean function associated with the
% switching dynamical SIM model.
% ARG sdsimKern : switching dynamical SIM kernel structure with the
% parameters.
% ARG t : input times for which the mean is to be computed.
% RETURN term : the value of $g_d$.
%
% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio A. Alvarez, 2015
    
term  = -sdsimKern.decay*exp(-sdsimKern.decay*t);

