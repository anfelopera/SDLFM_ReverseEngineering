function [X, y, XTest, yTest] = mapLoadDataLocal(dataset, seedVal)

% MAPLOADDATA Load a mapping model dataset (e.g. classification, regression).
% FORMAT
% DESC loads a data set for a mapping modelling problem
% (e.g. classification or regression).
% ARG dataset : the name of the data set to be loaded.
% RETURN X : the training input data loaded in.
% RETURN y : the training target data loaded in.
% RETURN XTest : the test input data loaded in. If no test set is
% available it is empty.
% RETURN yTest : a test target data.
%
% COPYRIGHT : Mauricio A. Alvarez, 2011

% DATASETS

if nargin < 2
    seedVal = 1e5;
end

RandStream('mt19937ar', 'Seed', seedVal);
XTest = [];
yTest = [];
baseDir = datasetsDirectory;

switch dataset
    case 'earthworm1'
        data = load([baseDir 'earthworms.mat']);
        labels1  = data.labels(:,1);
        % Put predefined classes in lower case
        for i=1:length(labels1)
            labels1{i} = lower(labels1{i});
        end
        % Look for predefined classes
        classes  = unique(labels1);
        K = length(classes);
        indexPerClass= cell(K,1);
        nPerClass = zeros(K,1);
        for i = 1: K
            indexPerClass{i} = find(strcmp(labels1, classes{i}));
            nPerClass(i) = length(indexPerClass{i});
        end
        chemicalFeatures = data.features(:, [1:6 8:13 15:20]);
        biodiversityFeatures = data.features(:, 22:30);
        latitude = data.features(:, end);
        featuresClasses = cell(K,1);
        for i=1:K
            K2 = length(indexPerClass{i});
            tempo = zeros(length(indexPerClass{i})*3, 10);
            cont = 0;
            cont2 =  0;
            cont3 = 0;
            for j=1:3  % three different depths
                tempo(cont3+1:cont3+K2, 1:6) = chemicalFeatures(indexPerClass{i}, cont+1:cont+6);
                cont = cont + 6;
                tempo(cont3+1:cont3+K2, 7:9) = biodiversityFeatures(indexPerClass{i}, cont2+1:cont2+3);
                tempo(cont3+1:cont3+K2, 10) = latitude(indexPerClass{i}, :);
                cont2 = cont2 + 3;
                cont3 = cont3 + K2;
            end
            featuresClasses{i} =  tempo;
        end
        % Further grouping according to the expert
        indexclass{1} = 2:8; % 'cafe caturra';'cafe con sombrio';'cafe/citricos'
        % 'cafe/platano';'cafe/platano/citrico';'cafe/platano/yuca'
        % 'cafetal variedad caturra asociado con yuca'
        indexclass{2} = 11:12; % 'guadua'; 'guadua/cacao'
        indexclass{3} = [10 13 16]; % 'frutales'; 'heliconia'; 'platano'
        indexclass{4} = [9 15]; % 'cana panelera'; 'pasto de corte'
        indexclass{5} = 1; % 'barbecho'
        indexclass{6} = 14; % 'pastizal'
        indexclass{7} = 17; % 'relicto de selva'
        indexclass{8} = 18; % 'variedad colombia'
        newK = 8;
        featuresClasses2 = featuresClasses;
        featuresClasses = cell(newK, 1);
        for i =1:newK
            featuresClasses{i} = cell2mat(featuresClasses2(indexclass{i}));
            
        end
        % Here, we only choose classes with N > 20.
        indexClassesFinalChosen = zeros(2,1);
        cont = 0;
        for i=1:newK
            if size(featuresClasses{i},1) > 20
                cont = cont + 1;
                indexClassesFinalChosen(cont) =  i;
            end
        end
        K2 = length(indexClassesFinalChosen);
        dataPerClass = zeros(K2,1);
        for i=1:K2
            dataPerClass(i) = size(featuresClasses{indexClassesFinalChosen(i)},1);
        end
        N = sum(dataPerClass);
        X = zeros(N, 10);
        y = zeros(N, 1);
        startOne = 1;
        endOne = 0;
        for i=1:K2
            endOne =  endOne + dataPerClass(i);
            X(startOne:endOne, :) = featuresClasses{indexClassesFinalChosen(i)};
            y(startOne:endOne, :) = i;
            startOne = endOne + 1;
        end
        XTest = classes(indexClassesFinalChosen);
    case 'earthworm2'
        data = load([baseDir 'earthworms.mat']);
        labels1  = data.labels(:,1);
        % Put predefined classes in lower case
        for i=1:length(labels1)
            labels1{i} = lower(labels1{i});
        end
        % Look for predefined classes
        classes  = unique(labels1);
        K = length(classes);
        indexPerClass= cell(K,1);
        nPerClass = zeros(K,1);
        for i = 1: K
            indexPerClass{i} = find(strcmp(labels1, classes{i}));
            nPerClass(i) = length(indexPerClass{i});
        end
        chemicalFeatures = data.features(:, [1:6 8:13 15:20]);
        biodiversityFeatures = data.features(:, 22:30);
        latitude = data.features(:, end);
        featuresClasses = cell(K,1);
        for i=1:K
            featuresClasses{i} = [chemicalFeatures(indexPerClass{i}, :) ...
                biodiversityFeatures(indexPerClass{i}, :) latitude(indexPerClass{i},:)];
        end
        % Further grouping according to the expert
        indexclass{1} = 2:8; % 'cafe caturra';'cafe con sombrio';'cafe/citricos'
        % 'cafe/platano';'cafe/platano/citrico';'cafe/platano/yuca'
        % 'cafetal variedad caturra asociado con yuca'
        indexclass{2} = 11:12; % 'guadua'; 'guadua/cacao'
        indexclass{3} = [10 13 16]; % 'frutales'; 'heliconia'; 'platano'
        indexclass{4} = [9 15]; % 'cana panelera'; 'pasto de corte'
        indexclass{5} = 1; % 'barbecho'
        indexclass{6} = 14; % 'pastizal'
        indexclass{7} = 17; % 'relicto de selva'
        indexclass{8} = 18; % 'variedad colombia'
        newK = 8;
        featuresClasses2 = featuresClasses;
        featuresClasses = cell(newK, 1);
        for i =1:newK
            featuresClasses{i} = cell2mat(featuresClasses2(indexclass{i}));
            
        end
        % Here, we only choose classes with N > 20.
        indexClassesFinalChosen = zeros(2,1);
        cont = 0;
        for i=1:newK
            if size(featuresClasses{i},1) > 10
                cont = cont + 1;
                indexClassesFinalChosen(cont) =  i;
            end
        end
        K2 = length(indexClassesFinalChosen);
        dataPerClass = zeros(K2,1);
        for i=1:K2
            dataPerClass(i) = size(featuresClasses{indexClassesFinalChosen(i)},1);
        end
        N = sum(dataPerClass);
        X = zeros(N, size(featuresClasses{1},2));
        y = zeros(N, 1);
        startOne = 1;
        endOne = 0;
        for i=1:K2
            endOne =  endOne + dataPerClass(i);
            X(startOne:endOne, :) = featuresClasses{indexClassesFinalChosen(i)};
            y(startOne:endOne, :) = i;
            startOne = endOne + 1;
        end
        XTest = classes(indexClassesFinalChosen);
        
    case 'earthwormsOldLady'
        data = load([baseDir 'earthwormsOldLady.mat']);
        labels1  = data.labels(:,1);
        % Put predefined classes in lower case
        for i=1:length(labels1)
            labels1{i} = lower(labels1{i});
        end
        % Look for predefined classes
        classes  = unique(labels1);
        K = length(classes);
        indexPerClass= cell(K,1);
        nPerClass = zeros(K,1);
        for i = 1: K
            indexPerClass{i} = find(strcmp(labels1, classes{i}));
            nPerClass(i) = length(indexPerClass{i});
        end
        chemicalFeatures = data.features(:, 1:24);
        biodiversityFeatures = data.features(:, 25:42);
        altitude = data.features(:, end);
        featuresClasses = cell(K,1);
        for i=1:K
            featuresClasses{i} = [chemicalFeatures(indexPerClass{i}, :) ...
                biodiversityFeatures(indexPerClass{i}, :) altitude(indexPerClass{i},:)];
        end
        % Further grouping according to the expert
        indexclass{1} = 1; % 'barbecho';
        indexclass{2} = 2:3; % 'cafetal arabigo'; 'cafetal asociado'
        indexclass{3} = 4; % 'cafetal variedad col'
        indexclass{4} = [5 7 11 13]; % 'cana'; 'pasto de corte'
        indexclass{5} = [6 8]; % 'eucalipto'; 'cipres'
        indexclass{6} = 9; % 'guadua'
        indexclass{7} = 10; % 'pastizal'
        indexclass{8} = 12; % 'relicto de selva'
        newK = 8;
        newClassesNames = cell(newK,1);
        for i = 1:newK
            concat = {};
            for j =1:length(indexclass{i})
                concat = {concat{:}, classes{indexclass{i}(j)}};
            end
            newClassesNames{i} = concat;
        end
        featuresClasses2 = featuresClasses;
        featuresClasses = cell(newK, 1);
        for i =1:newK
            featuresClasses{i} = cell2mat(featuresClasses2(indexclass{i}));
            
        end
        % Here, we only choose classes with N > 20.
        indexClassesFinalChosen = zeros(2,1);
        cont = 0;
        for i=1:newK
            if size(featuresClasses{i},1) > 20
                cont = cont + 1;
                indexClassesFinalChosen(cont) =  i;
            end
        end
        K2 = length(indexClassesFinalChosen);
        dataPerClass = zeros(K2,1);
        for i=1:K2
            dataPerClass(i) = size(featuresClasses{indexClassesFinalChosen(i)},1);
        end
        N = sum(dataPerClass);
        X = zeros(N, size(featuresClasses{1},2));
        y = zeros(N, 1);
        startOne = 1;
        endOne = 0;
        for i=1:K2
            endOne =  endOne + dataPerClass(i);
            X(startOne:endOne, :) = featuresClasses{indexClassesFinalChosen(i)};
            y(startOne:endOne, :) = i;
            startOne = endOne + 1;
        end
        XTest{1} = newClassesNames(indexClassesFinalChosen);
        XTest{2} = data.namesFeatures;
        
    case 'earthwormsFull'
        
        data = load([baseDir 'earthwormsFull.mat']);
        labels1  = data.labels(:,1);
        % Put predefined classes in lower case
        for i=1:length(labels1)
            labels1{i} = lower(labels1{i});
        end
        % Look for predefined classes
        classes  = unique(labels1);
        K = length(classes);
        indexPerClass= cell(K,1);
        nPerClass = zeros(K,1);
        for i = 1: K
            indexPerClass{i} = find(strcmp(labels1, classes{i}));
            nPerClass(i) = length(indexPerClass{i});
        end
        chemicalFeatures = data.features(:, 1:21);
        biodiversityFeatures = data.features(:, 22:27);
        altitude = data.features(:, end);
        featuresClasses = cell(K,1);
        for i=1:K
            featuresClasses{i} = [chemicalFeatures(indexPerClass{i}, :) ...
                biodiversityFeatures(indexPerClass{i}, :) altitude(indexPerClass{i},:)];
        end
        % Further grouping according to the expert
        indexclass{1} = 1; genclassNames{1} =  'barbecho'; % 'barbecho';
        indexclass{2} = 2:10; genclassNames{2} = 'cafe asociado';% cafe asociado: café con sombrío, café caturra,
        % café / cítricos, café / plátano, café / plátano / cítricos, café / plátano  / yuca,
        % cafetal var. caturra asociado Yuca.
        indexclass{3} = [11 25]; genclassNames{3} = 'cafetal variedad colombia';% 'cafetal variedad col'
        indexclass{4} = [12 13 15 17 20 22 23 26]; genclassNames{4} = 'cultivos temporales'; % cultivos temporales: cana, cana panelera,
        % citricos, frutales, heliconia, pasto de corte, platano, yuca
        indexclass{5} = [14 16]; genclassNames{5} = 'plantaciones';% plantaciones: 'eucalipto-lulo'; 'cipres'
        indexclass{6} = [18 19]; genclassNames{6} = 'guadua'; % guadua: 'guadua', 'guadua/cacao'
        indexclass{7} = 21; genclassNames{7} = 'pastizal'; % 'pastizal'
        indexclass{8} = 24; genclassNames{8} = 'relicto de selva'; % 'relicto de selva'
        newK = 8;
        newClassesNames = cell(newK,1);
        for i = 1:newK
            concat = {};
            for j =1:length(indexclass{i})
                concat = {concat{:}, classes{indexclass{i}(j)}};
            end
            newClassesNames{i} = concat;
        end
        featuresClasses2 = featuresClasses;
        featuresClasses = cell(newK, 1);
        for i =1:newK
            featuresClasses{i} = cell2mat(featuresClasses2(indexclass{i}));
            
        end
        % Here, we only choose classes with N > 20.
        indexClassesFinalChosen = zeros(2,1);
        cont = 0;
        for i=1:newK
            if size(featuresClasses{i},1) > 10
                cont = cont + 1;
                indexClassesFinalChosen(cont) =  i;
            end
        end
        K2 = length(indexClassesFinalChosen);
        dataPerClass = zeros(K2,1);
        for i=1:K2
            dataPerClass(i) = size(featuresClasses{indexClassesFinalChosen(i)},1);
        end
        N = sum(dataPerClass);
        X = zeros(N, size(featuresClasses{1},2));
        y = zeros(N, 1);
        startOne = 1;
        endOne = 0;
        for i=1:K2
            endOne =  endOne + dataPerClass(i);
            X(startOne:endOne, :) = featuresClasses{indexClassesFinalChosen(i)};
            y(startOne:endOne, :) = i;
            startOne = endOne + 1;
        end
        XTest{1} = genclassNames(indexClassesFinalChosen);
        XTest{2} = data.namesFeatures;
        
        
        
    case 'dbsUPVAdapWave1'
        data = load([baseDir 'dbsUPVAdapWave.mat']);
        features = data.Data_MTL{1}';
        labels = data.Data_MTL{2};
        %tasks = data.Data_MTL{3};
        % Look for predefined classes
        classes  = unique(labels);
        K = length(classes);
        indexPerClass= cell(K,1);
        nPerClass = zeros(K,1);
        for i = 1: K
            indexPerClass{i} = find( labels == classes(i));
            nPerClass(i) = length(indexPerClass{i});
        end
        N = sum(nPerClass);
        X = zeros(N, size(features,2));
        y = zeros(N, 1);
        startOne = 1;
        endOne = 0;
        for i=1:K
            endOne =  endOne + nPerClass(i);
            X(startOne:endOne, :) = features(indexPerClass{i}, :);
            y(startOne:endOne, :) = classes(i);
            startOne = endOne + 1;
        end
    case 'dbsUPVAdapWave2'
        
        data = load([baseDir 'dbsUPVAdapWave.mat']);
        tasks = data.Data_MTL{3};
        features = [data.Data_MTL{1}' tasks];
        labels = data.Data_MTL{2};
        % Look for predefined classes
        classes  = unique(labels);
        K = length(classes);
        indexPerClass= cell(K,1);
        nPerClass = zeros(K,1);
        for i = 1: K
            indexPerClass{i} = find( labels == classes(i));
            nPerClass(i) = length(indexPerClass{i});
        end
        N = sum(nPerClass);
        X = zeros(N, size(features,2));
        y = zeros(N, 1);
        startOne = 1;
        endOne = 0;
        for i=1:K
            endOne =  endOne + nPerClass(i);
            X(startOne:endOne, :) = features(indexPerClass{i}, :);
            y(startOne:endOne, :) = classes(i);
            startOne = endOne + 1;
        end
        
    case 'dbsUPVAdapWaveMultiTask1'
        
        data = load([baseDir 'dbsUPVAdapWave.mat']);
        tasks = data.Data_MTL{3};
        features = data.Data_MTL{1}';
        labels = data.Data_MTL{2};
        difTasks = unique(tasks);
        D = length(difTasks);
        indexPerClass= cell(D,1);
        for i = 1:D
            indexPerClass{i} = find( tasks  == difTasks(i));
        end
        X = cell(1, D);
        y = cell(1, D);
        for i=1:D
            X{i} = features(indexPerClass{i}, :);
            y{i} = labels(indexPerClass{i}, :);
        end
        
    case 'dbsUPVWaveMultiTask1'
        
        data = load([baseDir 'dbsUPVWave.mat']);
        tasks = data.Data_MTL{3};
        features = data.Data_MTL{1}';
        labels = data.Data_MTL{2};
        difTasks = unique(tasks);
        D = length(difTasks);
        indexPerClass= cell(D,1);
        for i = 1:D
            indexPerClass{i} = find( tasks  == difTasks(i));
        end
        X = cell(1, D);
        y = cell(1, D);
        for i=1:D
            X{i} = features(indexPerClass{i}, :);
            y{i} = labels(indexPerClass{i}, :);
        end
        
    case 'dbsUPVSpikesMultiTask1'
        
        data = load([baseDir 'dbsUPVSpikes.mat']);
        tasks = data.Data_MTL{3};
        features = data.Data_MTL{1}';
        labels = data.Data_MTL{2};
        difTasks = unique(tasks);
        D = length(difTasks);
        indexPerClass= cell(D,1);
        for i = 1:D
            indexPerClass{i} = find( tasks  == difTasks(i));
        end
        X = cell(1, D);
        y = cell(1, D);
        for i=1:D
            X{i} = features(indexPerClass{i}, :);
            y{i} = labels(indexPerClass{i}, :);
        end
        
    case 'dbsUTPAdapWave1'
        
        data = load([baseDir 'dbsUTPAdapWave.mat']);
        features = data.Data_MTL{1}';
        labels = data.Data_MTL{2};
        %tasks = data.Data_MTL{3};
        % Look for predefined classes
        classes  = unique(labels);
        K = length(classes);
        indexPerClass= cell(K,1);
        nPerClass = zeros(K,1);
        for i = 1: K
            indexPerClass{i} = find( labels == classes(i));
            nPerClass(i) = length(indexPerClass{i});
        end
        N = sum(nPerClass);
        X = zeros(N, size(features,2));
        y = zeros(N, 1);
        startOne = 1;
        endOne = 0;
        for i=1:K
            endOne =  endOne + nPerClass(i);
            X(startOne:endOne, :) = features(indexPerClass{i}, :);
            y(startOne:endOne, :) = classes(i);
            startOne = endOne + 1;
        end
        
    case 'dbsUTPAdapWave2'
        
        data = load([baseDir 'dbsUTPAdapWave.mat']);
        tasks = data.Data_MTL{3};
        features = [data.Data_MTL{1}' tasks];
        labels = data.Data_MTL{2};
        % Look for predefined classes
        classes  = unique(labels);
        K = length(classes);
        indexPerClass= cell(K,1);
        nPerClass = zeros(K,1);
        for i = 1: K
            indexPerClass{i} = find( labels == classes(i));
            nPerClass(i) = length(indexPerClass{i});
        end
        N = sum(nPerClass);
        X = zeros(N, size(features,2));
        y = zeros(N, 1);
        startOne = 1;
        endOne = 0;
        for i=1:K
            endOne =  endOne + nPerClass(i);
            X(startOne:endOne, :) = features(indexPerClass{i}, :);
            y(startOne:endOne, :) = classes(i);
            startOne = endOne + 1;
        end
        
    case 'dbsUTPAdapWaveMultiTask1'
        data = load([baseDir 'dbsUTPAdapWave.mat']);
        tasks = data.Data_MTL{3};
        features = data.Data_MTL{1}';
        labels = data.Data_MTL{2};
        difTasks = unique(tasks);
        D = length(difTasks);
        indexPerClass= cell(D,1);
        for i = 1:D
            indexPerClass{i} = find( tasks  == difTasks(i));
        end
        X = cell(1, D);
        y = cell(1, D);
        for i=1:D
            X{i} = features(indexPerClass{i}, :);
            y{i} = labels(indexPerClass{i}, :);
        end
        
    case 'dbsUTPWaveMultiTask1'
        
        data = load([baseDir 'dbsUTPWave.mat']);
        tasks = data.Data_MTL{3};
        features = data.Data_MTL{1}';
        labels = data.Data_MTL{2};
        difTasks = unique(tasks);
        D = length(difTasks);
        indexPerClass= cell(D,1);
        for i = 1:D
            indexPerClass{i} = find( tasks  == difTasks(i));
        end
        X = cell(1, D);
        y = cell(1, D);
        for i=1:D
            X{i} = features(indexPerClass{i}, :);
            y{i} = labels(indexPerClass{i}, :);
        end
        
    case 'dbsUTPSpikesMultiTask1'
        
        data = load([baseDir 'dbsUTPSpikes.mat']);
        tasks = data.Data_MTL{3};
        features = data.Data_MTL{1}';
        labels = data.Data_MTL{2};
        difTasks = unique(tasks);
        D = length(difTasks);
        indexPerClass= cell(D,1);
        for i = 1:D
            indexPerClass{i} = find( tasks  == difTasks(i));
        end
        X = cell(1, D);
        y = cell(1, D);
        for i=1:D
            X{i} = features(indexPerClass{i}, :);
            y{i} = labels(indexPerClass{i}, :);
        end
        
    case 'dbEmotionsCohnCanadeStaticReduced'
        
        data = load([baseDir 'emotionCohnCanadeReduced.mat']);
        %nemotions = size(data.emoAlign,2);
        nemotions = 2;
        nsamples = size(data.emoAlign(1,1).secuencia,2);
        nfeatures = size(data.emoAlign(1,1).secuencia(1).datos,1);
        X = zeros(nemotions*nsamples, nfeatures);
        y = zeros(nemotions*nsamples, 1);
        cont = 0;
        for i=1:nemotions         
            for j=1:nsamples
                cont = cont + 1;
                X(cont, :) = (data.emoAlign(i).secuencia(j).datos(:,end))';
                y(cont) = i;
            end
        end
        XTest = data.namesEmo(1:nemotions);
        
    case 'dbEmotionsCohnCanadeDynamicReduced'
        
        data = load([baseDir 'emotionCohnCanadeReduced.mat']);
        nemotions = size(data.emoAlign,2);
        nreps = size(data.emoAlign(1,1).secuencia,2);
        [noutputs, nlent] = size(data.emoAlign(1,1).secuencia(1).datos);
        X = cell(nemotions, 1);
        y = cell(nemotions, 1);
        for i=1:nemotions
            for j=1:nreps                
                X{i}{j} = mat2cell(repmat((1:nlent)', 1, noutputs), ...
                    nlent, ones(noutputs,1));                
                y{i}{j} = mat2cell((data.emoAlign(i).secuencia(j).datos)', ...
                    nlent, ones(noutputs,1));
            end
        end
        XTest = data.namesEmo;
        
    otherwise
        error('Database not recognized')
        
end
