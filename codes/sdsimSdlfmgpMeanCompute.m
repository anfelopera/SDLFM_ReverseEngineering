function m = sdsimSdlfmgpMeanCompute(meanFunction, X, varargin)

% SDSIMMEANCOMPUTE Give the output of the SDSIM mean function model for given X.
% FORMAT
% DESC gives the output of the sim mean function model for a given input X.
% ARG model : structure specifying the model.
% ARG X : input location(s) for which output is to be computed.
% RETURN Y : output location(s) corresponding to given input
% locations.
%
% SEEALSO : simMeanCreate
%
% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio A. Alvarez, 2015


% MULTIGP

nlf = varargin{1};
startVal=1;
endVal=0;

for i = nlf+1:nlf+length(X),
    endVal = endVal + size(X{i-nlf}, 1);    
    m(startVal:endVal, 1) = meanFunction.basal(i-nlf)/meanFunction.decay(i-nlf) * ...
        (1-exp(-meanFunction.decay(i-nlf)*X{i-nlf}));
    startVal = endVal+1;
end

% spVector = cumsum(meanFunction.switchingTimes);
% for i = nlf+1:length(X)+1,
%     mtemp = [];
%     endVal = endVal + size(X{i-nlf}, 1);    
%     for j=1:meanFunction.nIntervals-1
%         newt1 = X{i-nlf}(X{i-nlf}>=spVector(j) & X{i-nlf}<spVector(j+1));
%         mtemp = [mtemp 1-exp(-meanFunction.decay(i-nlf)*(newt1 - spVector(j)))'];
%     end
%     newt1 = X{i-nlf}(X{i-nlf}>=spVector(end));
%     mtemp = [mtemp 1-exp(-meanFunction.decay(i-nlf)*(newt1 - spVector(end)))'];
%     m(startVal:endVal, 1) = meanFunction.basal(i-nlf)/meanFunction.decay(i-nlf) *mtemp;
%     startVal = endVal+1;
% end