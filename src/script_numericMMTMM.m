
idx_solution = 1;
betaRange = linspace(0, 1, nBetaPoints);
alfaRange = linspace(0, 1, nAlphaPoints);
[m_alphaGrid, m_betaGrid] = meshgrid(alfaRange, betaRange);
for i = 1:size(inspectedPaths, 2)
    inspectedPaths(1, i).idx = i;
end

timeStartGlobal = tic;

% the main loop through frequencies
idx_frequency = linspace(1, 1001, nFrequencySteps);
idx_frequency = round(idx_frequency, 0);
for idx_currentFrequency = idx_frequency%1:frequencyStep:1001
    timeStartFrequencyLoop = tic;
    k0 = (2*pi*frequencies(idx_currentFrequency))/3e8;
    m_alphaUnnormalized = m_alphaGrid * k0;
    
    % sweep through multiple paths (each defined in init_script.m and added to vector in the main script)
    for s_path = inspectedPaths 
        m_betaUnnormalized = m_betaGrid * (pi / s_path.pathLength);
        
        % searching in complex plane
        for idx_beta = 1:nBetaPoints % sweep through BETA
            for idx_alpha = 1:nAlphaPoints% sweep throuth ALPHA
                k = m_betaUnnormalized(idx_beta, idx_alpha) - 1j*m_alphaUnnormalized(idx_beta, idx_alpha);
                
                % calculation of lambda matrix (eigenvalue matrix), used in
                % characteristic equation below
                m_lambda = getLambda(s_path, k, periodicity, nModes); 
                
                % calculating and storing results of haracteristic equation
                % for each frequency point, for each path
                s_path.m_solutionsPerOneFrequency(idx_beta, idx_alpha) = det(m_tParameters(:, :, idx_currentFrequency) - m_lambda);
            end
        end 
        % finding regional minimas in the currently solved path
        m_regionalMins = imregionalmin(db(abs(s_path.m_solutionsPerOneFrequency), 8));

        % Erasing numerical artifacts
        m_regionalMins(1:4, :) = 0;
        
        % transforming found minima into indexes
        idx_solutions = find(m_regionalMins);

        % normalization of beta and alpha values
        c_allSolutions{s_path.idx, idx_solution} = m_betaUnnormalized(idx_solutions)*(s_path.pathLength/pi) - 1j*m_alphaUnnormalized(idx_solutions)/k0;
        
    end
    % frequency is stored in the last line of cells
    c_allSolutions{nPaths+1, idx_solution} = frequencies(idx_currentFrequency);
    
    idx_solution = idx_solution + 1;
    clc
    fprintf("Paths: %d, beta steps: %d/%d; elapsed time per beta step: %.2f s\n", numel(inspectedPaths), idx_solution-1, nFrequencySteps, toc(timeStartFrequencyLoop))
end
fprintf("Done!\nTotal elapsed time: %.1f sec.\n", toc(timeStartGlobal))


clear idx_solution betaRange alfaRange idx_solutions m_lambda idx_beta column...
m_betaGrid m_alphaGrid m_regionalMins m_alphaUnnormalized m_betaUnnormalized