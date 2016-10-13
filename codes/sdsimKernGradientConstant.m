function [generalConstGrad, generalConst] = sdsimKernGradientConstant(nIntervals, ...
    simKern1, simKern2, spVector)

% SDSIMKERNGRADIENTCONSTANT Gradients for constants for the SDSIM kernel
% FORMAT
% DESC computes necessary gradients of the constants computed with
% sdsimKernComputeConstant.m
% ARG nIntervals : number of switching intervals in the kernel
% ARG simKern1 : structure containing the parameters of system 1
% ARG simkern2 : structure containing the parameters of system 2
% ARG spVector : vector containing the switching time values
% RETURN generalConstGrad : a cell containing the relative derivatives
% RETURN generalConstant : a cell containing the necessary constants for
% computing the kernel in the switching intervals.
%
% COPYRIGHT : Andres F. Lopez-Lopera and Mauricio A. Alvarez, 2015

% KERN

generalConst = cell(nIntervals);
generalGradDecay = cell(nIntervals);
generalGradSPoint = cell(nIntervals);

for i=1:nIntervals
    for j =1:i-1
        if i - j>=2
            c1 = sdsimMeanCompute(simKern1, spVector(i) - spVector(i-1));
            c2 = sdsimMeanCompute(simKern2, spVector(i) - spVector(i-1));
            g1 = sdsimvMeanCompute(simKern1, spVector(i) - spVector(i-1));
            g2 = sdsimvMeanCompute(simKern2, spVector(i) - spVector(i-1));
            % Derivatives wrt the decay so that we can get the
            % derivatives wrt paramters in a further step.
            gc1Decay = sdsimMeanGradient(simKern1, spVector(i) - spVector(i-1));
            gc2Decay = sdsimMeanGradient(simKern2, spVector(i) - spVector(i-1));
            gg1Decay = sdsimvMeanGradient(simKern1, spVector(i) - spVector(i-1));
            gg2Decay = sdsimvMeanGradient(simKern2, spVector(i) - spVector(i-1));
            % Derivatives wrt the switching points
            gc1sp = sdsimvMeanCompute(simKern1, spVector(i) - spVector(i-1));
            gc2sp = sdsimvMeanCompute(simKern2, spVector(i) - spVector(i-1));
            if i - j == 2
                % Computation of the constants
                generalConst{i,j}(1,1) = c1;
                generalConst{i,j}(2,1) = g1;
                generalConst{j,i}(1,1) = c2;
                generalConst{j,i}(2,1) = g2;
                % Gradient wrt the decay for simKern1
                generalGradDecay{i,j}(1,1) = gc1Decay;
                % Gradient wrt alpha for simKern2
                generalGradDecay{j,i}(1,1) = gc2Decay;
                % Gradient wrt switching points for simKern1
                generalGradSPoint{i,j}(1,1) = gc1sp; generalGradSPoint{i,j}(1,2) = -gc1sp;
                % Gradient wrt switching points for simKern2
                generalGradSPoint{j,i}(1,1) = gc2sp; generalGradSPoint{j,i}(1,2) = -gc2sp; % f1 wrt t_i and t_{i-1}
                
            else
                % Derivatives wrt to simKern1
                % Allocate values to make callings to
                % computeLocalDerivative short
                f1 = generalConst{i-1,j}(1,1);
                gf1Decay = generalGradDecay{i-1,j}(1,1);
                % Compute constants related to simKern1
                generalConst{i,j}(1,1) = c1*generalConst{i-1,j}(1,1);
                generalConst{i,j}(2,1) = g1*generalConst{i-1,j}(1,1);
                % Gradient wrt decay for simKern1
                generalGradDecay{i,j}(1,1) = computeLocalDerivative(c1, f1, gc1Decay, gf1Decay);
                generalGradDecay{i,j}(2,1) = computeLocalDerivative(g1, f1, gg1Decay, gf1Decay);
                % Gradient wrt to switching point t_i and t_{i-1}
                gf1sp = generalGradSPoint{i-1,j}(1,1);
                % First the derivative wrt to t_i
                generalGradSPoint{i,j}(1,1) = gc1sp*f1;
                % Second the derivative wrt to t_{i-1}
                generalGradSPoint{i,j}(1,2) = computeLocalDerivative(c1, f1, -gc1sp, gf1sp);
                % Now compute all other derivatives wrt to the last
                % switching points
                maxL = i - j;
                cont = 1;
                for k = 3:maxL,
                    cont = cont + 1;
                    generalGradSPoint{i,j}(1,k) = c1*generalGradSPoint{i-1,j}(1,cont);
                end
                % Derivatives wrt to simKern2
                % Allocate values to make callings to
                % computeLocalDerivative short
                f1 = generalConst{j,i-1}(1,1);
                gf1Decay = generalGradDecay{j,i-1}(1,1);
                % Compute constants related to simKern2
                generalConst{j,i}(1,1) = c2*generalConst{j,i-1}(1,1);
                generalConst{j,i}(2,1) = g2*generalConst{j,i-1}(1,1);
                % Gradient wrt decay for simKern2
                generalGradDecay{j,i}(1,1) = computeLocalDerivative(c2, f1, gc2Decay, gf1Decay);
                % Gradient wrt to switching point t_i and t_{i-1}
                gf1sp = generalGradSPoint{j, i-1}(1,1);
                % First the derivative wrt to t_i
                generalGradSPoint{j,i}(1,1) = gc2sp*f1;
                % Second the derivative wrt to t_{i-1}
                generalGradSPoint{j,i}(1,2) = computeLocalDerivative(c2, f1, -gc2sp, gf1sp);
                % Now compute all other derivatives wrt to the last
                % switching points
                maxL = i - j;
                cont = 1;
                for k = 3:maxL,
                    cont = cont + 1;
                    generalGradSPoint{j,i}(1,k) = c2*generalGradSPoint{j,i-1}(1,cont);
                    generalGradSPoint{j,i}(2,k) = g2*generalGradSPoint{j,i-1}(1,cont);
                end
                
                
            end
        end
    end
end

generalConstGrad{1} = generalGradDecay;
generalConstGrad{2} = generalGradSPoint;

function gradtot = computeLocalDerivative(a, b, grada, gradb)

gradtot = a*gradb + grada*b;
