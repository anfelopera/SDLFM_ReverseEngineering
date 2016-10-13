function k = sdmultisdsimKernDiagCompute(kern, x)

% SDMULTISDSIMKERNDIAGCOMPUTE Compute diagonal of SDMULTI for SDSIM kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the switching
% dynamical multiple output block kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : sdmultiKernParamInit, kernDiagCompute, sdmultiKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Pei Gao, 2007

% MODIFICATIONS : Andres F. Lopez Lopera and Mauricio A. Alvarez, 2015

% SDSIMGP

% It actually calls the same function that already exists for multiKern

% Divide the matrix of initial conditions in a matrix of proper dimensions
% for positions and velocities

  kyyTemp = kern.KIC(1:kern.numPositions, 1:kern.numPositions);

  % We repeat the initial condition matrices, to make easier the process of
  % passing them to the sdmultiKernBlock function
  howmany = 1 + kern.includeVel + kern.includeAccel;
  kyy = repmat(kyyTemp, howmany, howmany);
    
  if iscell(x)
    dim = 0;
    for i = 1:length(x)
      dim = dim + size(x{i},1);
    end
    k = zeros(dim, 1);
    startVal = 1;
    endVal = size(x{1},1);
    for i = 1:length(kern.comp)
      covIC(1,1) = kyy(i,i);     
      if ~isempty(x{i})
        k(startVal:endVal) = sdkernDiagCompute(kern.comp{i}, x{i}, covIC);
      end
      startVal = endVal + 1;
      if i+1 <= length(kern.comp)
        endVal = endVal + size(x{i+1},1);
      end
    end
  else
    k = zeros(size(x, 1)*kern.numBlocks, 1);
    endVal = size(x, 1);
    startVal = 1;
    for i = 1:length(kern.comp)
      covIC(1,1) = kyy(i,i);
      k(startVal:endVal, 1)  = sdkernDiagCompute(kern.comp{i}, x, covIC);
      startVal = endVal + 1;
      endVal = endVal + size(x, 1);
    end
  end
  