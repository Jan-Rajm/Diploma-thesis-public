%% The main script of MMTMM (MultiModal Transfer Matrix Method)
% Convention: - prefix s_ represents structures
%             - prefix o_ represents objects
%             - prefix m_ represents matrices
%             - prefix c_ represents cells
%             - prefix script_ represents calling of another script
%             - prefix idx_ represents vectors of logic values as indices
%             - no prefix is used for general single values and general vectors
%             - single-letter variables are used only for variables, related
%               to physical quantities, such as wavenumber : k
% CammelCase is used for variable names
% Regarding names of paths, "G" shall be understood as "Gamma" (GX = Gamma-X)
clc
close all
clear all
%% User s_settings

nFrequencySteps = 100;
nBetaPoints = 150;
nAlphaPoints = 30;

simulationProfile = "3D_cube"; %"2D_square" "3D_cube" "2D_hexagonal" "lens"

% setting structures, unit cell related parameters and paths are defined inside
script_initialization; 

% Uncomment to overwrite paths defined in init_script
% inspectedPaths = [Gx Xm Mg1];% Xm Mg1];
% inspectedPaths = [Gx];
% Setting structure, determines which processes will be executed
s_settings = s_numericMMTMM;

% Postprocessing and graphical output
addEigenmode = "";% name of .txt file, containing eigenmode plot
resultsToMerge = ""; % create dispersion diagram from multiple stored results
%% Process of MMTMM
nPaths = size(inspectedPaths, 2);
if ~s_settings.processSaved
    % definition of simulation variables
    frequencyStep = 1001 / nFrequencySteps; % 1001 points set in CST while exporting touchstone files
    s_simulationData = sparameters(touchstoneFile);
    nPorts = s_simulationData.NumPorts;
    nModes = nPorts / (2 * numel(periodicity));  
    frequencies = s_simulationData.Frequencies;
    
    % preallocation of data matrices, definition of matrix variables 
    m_sParameters = s_simulationData.Parameters;
    m_tParameters = s2abcd(m_sParameters, s_simulationData.Impedance);
    c_allSolutions = cell(nPaths + 1, nFrequencySteps);
    
    % scripts of MMTMM itself
    if s_settings.numeric
        script_numericMMTMM;
    end
    
    if s_settings.linearized
        script_linearizedMMTMM;
    end
    
    if s_settings.save
        script_saveResults;
    end
end

%% Postprocessing
if s_settings.process_saved
    c_allSolutions = [];
    for resultFile = resultsToMerge
        storedResults = open(resultFile).c_allSolutions;
        c_allSolutions = [c_allSolutions storedResults];
    end
end

% defined as cells, because modes can differ in count of points
% rows = path (e.g. GX)
c_resultBeta = {}; % all beta values
c_resultAlpha = {}; % all alpha values
c_resultFrequency = {}; % all frequencies

% In this script are obtained data rearranged into more convenient format
script_processPlotData;


% Shift different paths to create overall view of dispersion diagram
for i = 1:nPaths
    c_resultBeta{i, :} = c_resultBeta{i, :} + i - 1;
end

% If an eigenmode plot data file is defined, it will be added to dispersion
% diagram. Following script reads ASCII plot data file and rearranges into
% more convenient format.
if addEigenmode ~= ""
    script_processEigenmode;
end


%% Ploting data
% Beta plot
figure;
hold on
if addEigenmode ~= ""
    % Due to the way of adding multiple plots into one figure, this
    % workaround is necessary to display legend correctly
    plot(0, NaN, "bx")
    plot(0, NaN, "ro")
end
for path = 1:nPaths
    plot(c_resultBeta{path, :}, c_resultFrequency{path, :}*1e-9, "bx")
end

grid on

% Configuration of x-axis in order to display critical points instead of
% numbers
o_xAxisProperties = get(gca, 'XAxis');
o_xAxisProperties.TickLabelInterpreter = 'latex';
if isempty(globalXticks)
    xticks(0:nPaths+1)
else
    xticks(globalXticks)
end
xticklabels(dispersionDiagramXLabels)
if exist("yTicks", "var")
    yticks(yTicks)
else
    yticks(0:5:50)
end
hold off
xlabel("Phase shift $\beta p/\pi$", "Interpreter", "latex")
ylabel("Frequency [GHz]", "Interpreter", "latex")
ylim(yLimitsBeta)

if addEigenmode ~= ""
    hold on
    for i = 1:size(eigenvalueX, 1)
        normalizedEigenvalX = eigenvalueX(i, :);
        normalizedEigenvalX = normalizedEigenvalX / max(normalizedEigenvalX);
        normalizedEigenvalX = normalizedEigenvalX * nPaths;

        plot(normalizedEigenvalX, eigenvalueY(i, :), "ro")
    legend("MMTMM", "Eigenmode")
    end
    hold off
end

% Alpha plot
figure
hold on
for path = 1:nPaths
    Alpha_plt = plot(c_resultAlpha{path, :}, c_resultFrequency{path, :}*1e-9, "x");
end
hold off
xlabel("Attenuation $\alpha/k_{0}$", "Interpreter", "latex")
ylabel("Frequency [GHz]", "Interpreter", "latex")
xlim(xLimitsAlpha)
ylim(yLimitsAlpha)
grid on

for idx_path = 1:nPaths
    currentPath = inspectedPaths(1, idx_path);
    tempFrequency = c_resultFrequency{idx_path, 1}; %2 pi f / c
    tempBeta = c_resultBeta{idx_path, 1};
    k0 = 2 * pi * tempFrequency / 3e8;

    % Unnormalization of beta values
    tempBeta = tempBeta ./ max(tempBeta);
    tempBeta = tempBeta * (pi / currentPath.pathLength);
    
    n = tempBeta ./ k0;
    
    % Calculation of appropriate maximal plot vaulue
    yMax = max(n) / 0.5;
    yMax = round(yMax) * 0.5 + 0.5;
   
    figure
    plot(tempFrequency, n, "bx")
    title(currentPath.name)
    xlabel("Frequency [GHz]")
    ylabel("Refractive index [-]")
    ylim([0 yMax])
    grid on

end
