%% Demonstration script of the Periodicity of modes base method
%  Demonstration is conducted on eigenmode plot data, obtained for one path
%
%  Following script can unfold the first mode of given path for arbitrary
%  enlargement factor

clc
clear all
close all

enlargementFactor = 4;
supercellEigenmodeData = readmatrix("2d_" + num2str(enlargementFactor) + ".txt")';
referenceEigenmodeData = readmatrix("primitiveCell.txt")';
beta = supercellEigenmodeData(1, :);
frequency = supercellEigenmodeData(2, :);

% Hold lower boundary of beta value of the currently processed region.
% Upper boudaty beta value is assumed to be betaRegion + (180 / enlargementFactor)
betaRegion = 0;

% Hold lower boundary of frequency value of the currently processed region.
% Upper boudaty frequency value is assumed to be
% betaRegion + (modeMaxFrequency / enlargementFactor)
frequencyRegion = 0;

% Needs to be set manualy
modeMaxFrequency = 24.26;

figure
hold on
colorSupercell = "bx";
colorEstimatedPrimitiveCell = "r-o";
for i = 1:enlargementFactor
    idx_betaRegion = (beta >= betaRegion) & (beta <= betaRegion + 180);
    idx_frequencyRegion = (frequency >= frequencyRegion) & (frequency <= frequencyRegion + modeMaxFrequency/enlargementFactor);
    plot(beta(idx_betaRegion), frequency(idx_betaRegion), colorSupercell)
    idx = idx_betaRegion & idx_frequencyRegion;
    plot(beta(idx), frequency(idx), colorEstimatedPrimitiveCell)
    betaRegion = betaRegion + 180;
    frequencyRegion = frequencyRegion + modeMaxFrequency / enlargementFactor;
end
plot(referenceEigenmodeData(1, :) * enlargementFactor, referenceEigenmodeData(2, :), "g--")
hold off
grid on
xlabel("Phase shift $\beta p/\pi$", "Interpreter", "latex")
ylabel("Frequency [GHz]", "Interpreter", "latex")
