%% Data processing mode

% Predefined setting structures
% .numeric      : should be numeric variant of MMTMM be launched?
% .save         : should be obtained values of beta, alpha and frequency be stored?
% .processSaved : should stored results be loaded?

s_mergeResults.numeric = false;
s_mergeResults.linearized = false;
s_mergeResults.save = false;
s_mergeResults.processSaved = true;

s_numericMMTMM.numeric = true;
s_numericMMTMM.linearized = false;
s_numericMMTMM.save = false;
s_numericMMTMM.processSaved = false;

% Predefined matrix of solutions per one frequency, stores results of the
% characteristic equation
m_solution = NaN * ones(nBetaPoints, nAlphaPoints);

% Numerical artifacts may appear, when alpha reach too heigh values
% If this occures, both alpha value and corresponding beta value are
% excluded
alphaValuesFilter = 1e3;
%% Definitions of directions
%% 2D conventional structure, dielectric
if simulationProfile == "2D_square"
    touchstoneFile = "";                     % Define input touchstone file
    p = 0;                                   % Define periodicity in meters 
    periodicity = [p p];
    dispersionDias_gramXLabels = ["$\Gamma$", "X", "M", "$\Gamma$"];
    globalXticks = 0:numel(dispersionDias_gramXLabels) - 1;

    yLimitsBeta = [0 50];                   % Define limits of x and y axis
    xLimitsAlpha = [0 1];
    yLimitsAlpha = yLimitsBeta;

    s_s_gx.xScale = 0;
    s_s_gx.yScale = 1;
    s_gx.xOffset = 0;
    s_gx.yOffset = 0;
    s_gx.m_solutionsPerOneFrequency = m_solution;
    s_gx.name = "GX";
    s_gx.pathLength = p;
    s_gx.idx = -1;
    
    s_xm.xScale = 1;
    s_xm.yScale = 0;
    s_xm.xOffset = 0;
    s_xm.yOffset = pi/p;
    s_xm.m_solutionsPerOneFrequency = m_solution;
    s_xm.name = "XM";
    s_xm.pathLength = p;
    s_xm.idx = -1;

    s_mg.xScale = 1;
    s_mg.yScale = 1;
    s_mg.xOffset = -pi/p;
    s_mg.yOffset = -pi/p;
    s_mg.m_solutionsPerOneFrequency = m_solution;
    s_mg.name = "GM";
    s_mg.pathLength = p;
    s_mg.idx = -1;

    inspectedPaths = [s_gx s_xm s_mg];
end
%% HEX 2
% 2D structure, metalic + dielectric
if simulationProfile == "2D_hexagonal"
    touchstoneFile = "";
    dispersionDias_gramXLabels = ["$\Gamma$", "M", "K", "$\Gamma$"];
    px = 0e-3;
    py = 0e-3;
    periodicity = [px py];

    yLimitsBeta = [0 50];
    xLimitsAlpha = [0 1];
    yLimitsAlpha = yLimitsBeta;

    globalXticks = 0:numel(dispersionDias_gramXLabels) - 1;
    
    s_mg.xScale = 2/3;
    s_mg.yScale = 1;
    s_mg.xOffset = -1*(2/3)*(pi/px);
    s_mg.yOffset = -1*(2*pi/py);
    s_mg.m_solutionsPerOneFrequency = m_solution;
    s_mg.name = "MG";
    s_mg.pathLength = (0.2)*1.3*getHypotenuse(px*(2/3), 2*py);
    s_mg.idx = -1;

    s_gx.xScale = 0;
    s_gx.yScale = 2;
    s_gx.xOffset = 0;
    s_gx.yOffset = 0;
    s_gx.m_solutionsPerOneFrequency = m_solution;
    s_gx.name = "GX";
    s_gx.pathLength = py;
    s_gx.idx = -1;
    
    s_xm.xScale = 2/3;
    s_xm.yScale = 0;
    s_xm.xOffset = 0;
    s_xm.yOffset = 2*pi/py;
    s_xm.m_solutionsPerOneFrequency = m_solution;
    s_xm.name = "XM";
    s_xm.pathLength = px;
    s_xm.idx = -1;

    inspectedPaths = [s_gx s_xm s_mg];
end
%% 3D
% 3D structure, metalic
if simulationProfile == "3D_cube"
    touchstoneFile = "";
    dispersionDias_gramXLabels = ["$\Gamma$" "X" "M" "$\Gamma$" "R" "X,M" "R"];
    globalXticks = 0:numel(dispersionDias_gramXLabels) - 1;
    p = 0e-3; 
    periodicity = [p p p];
    
    alphaValuesFilter = 1e-3;

    yLimitsBeta = [0 50];
    xLimitsAlpha = [0 1];
    yLimitsAlpha = yLimitsBeta;

    s_gx.xScale = 1;
    s_gx.yScale = 0;
    s_gx.zScale = 0;
    s_gx.xOffset = 0;
    s_gx.yOffset = 0;
    s_gx.zOffset = 0;
    s_gx.m_solutionsPerOneFrequency = m_solution;
    s_gx.name = "GX";
    s_gx.pathLength = p;
    s_gx.idx = -1;

    s_xm.xScale = 0;
    s_xm.yScale = 1;
    s_xm.zScale = 0;
    s_xm.xOffset = pi/p;
    s_xm.yOffset = 0;
    s_xm.zOffset = 0;
    s_xm.m_solutionsPerOneFrequency = m_solution;
    s_xm.name = "XM";
    s_xm.pathLength = p;
    s_xm.idx = -1;
    
    s_mg.xScale = 1;
    s_mg.yScale = 1;
    s_mg.zScale = 0;
    s_mg.xOffset = -pi/p;
    s_mg.yOffset = -pi/p;
    s_mg.zOffset = 0;
    s_mg.m_solutionsPerOneFrequency = m_solution;
    s_mg.name = "MG";
    s_mg.pathLength = p;
    s_mg.idx = -1;
    
    s_gr.xScale = 1;
    s_gr.yScale = 1;
    s_gr.zScale = 1;
    s_gr.xOffset = 0;
    s_gr.yOffset = 0;
    s_gr.zOffset = 0;
    s_gr.m_solutionsPerOneFrequency = m_solution;
    s_gr.name = "GR";
    s_gr.pathLength = p;
    s_gr.idx = -1;
    
    s_rx.xScale = 0;
    s_rx.yScale = 1;
    s_rx.zScale = 1;
    s_rx.xOffset = pi/p;
    s_rx.yOffset = -pi/p;
    s_rx.zOffset = -pi/p;
    s_rx.m_solutionsPerOneFrequency = m_solution;
    s_rx.name = "RX";
    s_rx.pathLength = p;
    s_rx.idx = -1;

    s_mr.xScale = 0;
    s_mr.yScale = 0;
    s_mr.zScale = 1;
    s_mr.xOffset = pi/p;
    s_mr.yOffset = pi/p;
    s_mr.zOffset = 0;
    s_mr.m_solutionsPerOneFrequency = m_solution;
    s_mr.name = "MR";
    s_mr.pathLength = p;
    s_mr.idx = -1;

    inspectedPaths = [s_gx s_xm s_mg s_gr s_rx s_mr];
end