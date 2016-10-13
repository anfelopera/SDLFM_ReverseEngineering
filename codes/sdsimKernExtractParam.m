function [params, names] = sdsimKernExtractParam(kern)

% SDSIMKERNEXTRACTPARAM Extract parameters from the SDSIM kernel structure.
% FORMAT
% DESC Extract parameters from the switching dynamical SIM kernel structure
% into a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. 
%
% FORMAT
% DESC Extract parameters and their names from the switching dynamical
% SIM kernel structure.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. 
% RETURN names : cell array of strings containing parameter names.
%
% SEEALSO sdsimKernParamInit, sdsimKernExpandParam, kernExtractParam, 
%
% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio A. Alvarez, 2015

% KERN

params = [kern.decay, kern.inverseWidth(:)', kern.switchingTimes, kern.sensitivity(:)'];
if isfield(kern.options, 'priorS') && kern.options.priorS
    params = [params kern.priorSvariance(:)'];
end

if nargout > 1
    names = {'decay'};
    namesInvWidth = cell(kern.nIntervals*kern.nlfPerInt,1);
    namesSensitivities = cell(kern.nIntervals*kern.nlfPerInt,1);
    if isfield(kern.options, 'priorS') && kern.options.priorS
        namesPriorSvariances = cell(kern.nIntervals*kern.nlfPerInt,1);
    end
    if kern.nlfPerInt == 1
        for i=1:kern.nIntervals
            namesInvWidth{i} = ['inverse width interval ' num2str(i) '.'];
            namesSensitivities{i} = ['sensitivity interval ' num2str(i) '.'];
            if isfield(kern.options, 'priorS') && kern.options.priorS
                namesPriorSvariances{i} = ['priorSvariance interval ' num2str(i) '.'];
            end
        end
    else
        cont = 0;
        for i=1:kern.nlfPerInt
            for j=1:kern.nIntervals
                cont = cont + 1;
                namesInvWidth{cont} = ['inverse width ' num2str(i) '.' ' interval ' num2str(j) '.'];
                namesSensitivities{cont} = ['sensitivity ' num2str(i) '.' ' interval ' num2str(j) '.'];
                if isfield(kern.options, 'priorS') && kern.options.priorS
                    namesPriorSvariances{cont} = ['priorSvariance ' num2str(i) '.' ' interval ' num2str(j) '.'];
                end                
            end
        end
    end
    namesStimes = cell(kern.nIntervals, 1);
    for i=1:kern.nIntervals
        namesStimes{i} = ['switching point interval ' num2str(i) '.'];
    end
    names = {names{:}, namesInvWidth{:}, namesStimes{:}, ...
        namesSensitivities{:}};
    if isfield(kern.options, 'priorS') && kern.options.priorS
        names = {names{:} namesPriorSvariances{:}};
    end
end
