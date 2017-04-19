function  ll = multigpLogLikelihood( model)

% MULTIGPLOGLIKELIHOOD Compute the log likelihood of a MULTIGP.

% COPYRIGHT : Mauricio A Alvarez, 2008

% MODIFICATIONS : Andres F. Lopez-Lopera and Mauricio A. Alvarez, 2015

% MULTIGP

switch model.approx
    case 'ftc'
        dim = size(model.m, 1);
        ll = -dim*log(2*pi) - model.logDetK - model.m'*model.invK*model.m;
        ll = ll*0.5;
        if isfield(model.kern.comp{1, 1}.comp{1, 1}, 'priorS') ...
                && model.kern.comp{1, 1}.comp{1, 1}.priorS ;
            nS = model.nout*model.nlf*model.nIntervals;
            S  = []; 
            Svar = [];
            for k = 1:length(model.kern.comp{1, 1}.comp)
                S  = [S; model.kern.comp{1, 1}.comp{1, k}.sensitivity(:)];
                Svar = [Svar; model.kern.comp{1, 1}.comp{1, k}.priorSvariance(:)];
            end
            lls= -nS*log(2*pi) - sum(log(Svar)) - sum( (S.^2)./Svar );            
            ll = ll + 0.5*lls;
        end
        % MAP optimization for gpsim with more than two latent functions
        %         if strcmp(model.kernType, 'sim') && (model.nlf>1) && ...
        %                 strcmp(model.inference, 'map')
        %             ll = ll - sum(log(model.S));
        %         end
    case {'sor','dtc','fitc','pitc','dtcvar'}
        ll = spmultigpLogLikelihood( model);
end

