function tieInd = sdsimSdlfmgpTieParam(model, options)

% SDLFMSDLFMGPTIEPARAM Tie parameters for a sdlfmgp model with sdsim kernel
% FORMAT
% DESC Tie the parameters for a sdlfmgp model that uses a sdsim kernel.
% RETURN tieInd : cell with elements containing the indexes of parameters
% to tie.
% ARG model : model created
% ARG options : options for the model.
%
% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio A. Alvarez 2015

% SDLFMGP

tieInd = cell(model.nIntervals*(model.nlfPerInt+1),1);
% Tie the parameters for the inverse widths
cont = 0;
for i=1:model.nIntervals
    if model.nlfPerInt == 1
        indx = paramNameRegularExpressionLookup(model, ...
            ['inverse width interval ' num2str(i) '\.']); 
        cont = cont + 1;
        tieInd{cont} = indx;
    else
        for j=1:model.nlfPerInt
            indx = paramNameRegularExpressionLookup(model, ['.* inverse width ' num2str(j) '\.' ...
                ' interval ' num2str(i) '\.']); 
            cont = cont + 1;
            tieInd{cont} = indx;
        end
    end
end
% Tie the parameters of the switching times
for i=1:model.nIntervals
    indx = paramNameRegularExpressionLookup(model, ...
        ['.* switching point interval ' num2str(i) '\.']);  
    cont = cont + 1;
    tieInd{cont} = indx;   
end

if options.includeVel || options.includeAccel
    decayIndx = paramNameRegularExpressionLookup(model, '.* decay');
    for i=1:model.numPositions
        if options.includeVel && ~options.includeAccel
            indxdecay = [decayIndx(i) decayIndx(i+model.numPositions)];
        else
            indxdecay = [decayIndx(i) decayIndx(i+model.nout) decayIndx(i+2*model.numPositions)];
        end
        cont = cont + 1;
        tieInd{cont} = indxdecay;
    end
end

if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
    for i = 1:model.nout 
        indx = paramNameRegularExpressionLookup(model, [num2str(i) ' decay']);
        cont = cont + 1;
        tieInd{cont} = indx;        
    end
%     for i=1:model.nIntervals
%         indx = paramNameRegularExpressionLookup(model, ...
%                 ['.* switching point interval ' num2str(i) '\.']);  
%         cont = cont + 1;
%         tieInd{cont} = indx;
%    end
end
