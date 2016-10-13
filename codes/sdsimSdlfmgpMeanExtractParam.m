function [params, names] = sdsimSdlfmgpMeanExtractParam(meanFunction)

% SDSIMMEANEXTRACTPARAM Extract parameters from the SDSIM MEAN FUNCTION structure.
% FORMAT
% DESC Extract parameters from the mean funtion structure of the sdsim model
% into a vector of parameters for optimisation.
% ARG meanFunction : the mean function structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. 
%
% SEEALSO  sdsimSdlfmgpMeanCreate, sdsimMeanExpandParam, sdlfmgpKernCreate,
% sdsimkernExtractParam
%
% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio A. Alvarez, 2015

% MULTIGP

params = [meanFunction.basal' meanFunction.decay' meanFunction.switchingTimes];
if nargout > 1
    names = cell(1, 2*meanFunction.nout);
    for i=1:meanFunction.nout
        names{i} = ['sdsim ' num2str(i) ' basal'];
    end    
    for i=meanFunction.nout+1:2*meanFunction.nout        
        names{i} = ['sdsim ' num2str(i-meanFunction.nout) ' decay'];
    end 
    namesStimes = cell(meanFunction.nIntervals, 1);
    for i=1:meanFunction.nIntervals
        namesStimes{i} = ['switching point interval ' num2str(i) '.'];
    end
    names = {names{:}, namesStimes{:}};
end
