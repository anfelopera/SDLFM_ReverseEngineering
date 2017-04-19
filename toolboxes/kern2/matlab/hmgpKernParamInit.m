function kern = hmgpKernParamInit(kern)

% TO DO

if isfield(kern, 'options') && isfield(kern.options, 'nlf')
    kern.nlf = kern.options.nlf;
else
    kern.nlf = 1;
end

kernType = {'parametric', kern.options, kern.options.kernType};

if kern.nlf > 1
    kernType1 = cell(1, kern.nlf);
    for i = 1:kern.nlf
        kernType1{i} = kernType; 
    end
    kglobal = kernCreate(1, {'cmpnd', kernType1{:}});
else
    kglobal = kernCreate(1, kernType);
end









if isfield(kern, 'options') && isfield(kern.options, 'basicKernelType')
    kern.basicKernelType = kern.options.basicKernelType;
else
    kern.basicKernelType = 'gaussian';
end

if isfield(kern, 'options') && isfield(kern.options, 'rankCorregMatrix')
   kern.rankCorregMatrix = kern.options.rankCorregMatrix; 
else
   kern.rankCorregMatrix = 1; 
end

fhandle = str2func([kern.basicKernelType 'KernParamInit']);
if isfield(kern, 'options') && isfield(kern.options, 'isArd')
   kern = fhandle(kern, kern.options.isArd); 
else
   kern = fhandle(kern);
end

if isfield(kern, 'options') && isfield(kern.options, 'nout')
    kern.nout = kern.options.nout;
    kern.A = rand(kern.nout, kern.rankCorregMatrix); 
else
    % If the number of outputs is not provided, we assume is one
    kern.nout = 1;
    kern.A = rand(1, kern.rankCorregMatrix);
end
if kern.nout<kern.rankCorregMatrix
    error('The rank of the matrix is greater than the number of outputs')
end
kern.B = kern.A*kern.A';
% Count again the number of parameters
kern.nParamsBK = kern.nParams;
kern.nParams = kern.nParams + kern.rankCorregMatrix*kern.nout;
