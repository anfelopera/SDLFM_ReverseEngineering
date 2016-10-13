function g = sdsimSdlfmgpMeanGradient(meanFunction, varargin)

% SDSIMMEANGRADIENT Gradient of the parameters of the mean function in the
% multigp model with SDSIM kernel
% FORMAT
% DESC gives the gradient of the objective function for the parameters of
% the mean function in the multigp model with SDSIM kernel (first order
% differential equation).
% ARG meanFunction : mean function structure to optimise.
% ARG P1, P2, P3 ... : optional additional arguments.
% RETURN g : the gradient of the error function to be minimised.
% 
% SEEALSO : simMeanCreate, simMeanOut
%
% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio Alvarez, 2015

% MULTIGP

gmu = varargin{1}';
X   = meanFunction.X;
d   = length(X);
switchingTimes = meanFunction.switchingTimes;
q   = length(switchingTimes);

gB = zeros(d,1);
gD = zeros(d,1);
gsp= zeros(q,d);

spVector = cumsum(switchingTimes);

for k = 1:d
    cd  = exp(-meanFunction.decay(k)*X{k});
    gB(k) = gmu{k}*(1-cd)/meanFunction.decay(k);
    gD(k) = gmu{k}*(X{k}.*cd) - gB(k);
    gD(k) = gD(k)*meanFunction.basal(k)/meanFunction.decay(k);
    for j = 1:q
        gsp(k,j) = gmu{k}*(exp(-meanFunction.decay(k)*(X{k} - spVector(j))))*meanFunction.basal(k);        
    end
end
gSp = sum(gsp,1);
gSp = zeros(1,q);
g = [gB' gD' gSp];
