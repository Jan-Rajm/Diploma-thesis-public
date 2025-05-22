%% Merging solutions (= columns 1:end-1) of multiple paths into distiguished
%  cell arrays
for idx_path = 1:nPaths
    betaTemp = [];
    alphaTemp = [];
    tempFrequency = [];
    for solution = c_allSolutions % sweep through columns
        betaAppended  = abs(real(solution{idx_path, 1}'));
        alphaAppended = abs(imag(solution{idx_path, 1}'));

        idx_temp = alphaAppended > alphaValuesFilter;
        alphaAppended(idx_temp) = NaN;
        betaAppended(idx_temp) = NaN;

        betaTemp = [betaTemp betaAppended];
        alphaTemp = [alphaTemp alphaAppended];

        % both alphaAppended and betaAppended can store multiple values
        % thus count of frequency points must be 
        tempFrequency = [tempFrequency solution{end, 1} * ones(1, size(betaAppended, 2))];
    end
    c_resultBeta{end+1, 1} = betaTemp;
    c_resultAlpha{end+1, 1} = alphaTemp;
    c_resultFrequency{end+1, 1} = tempFrequency;
end

% clearing temporal variables, which will not be needed anymore
clear solution alphaTemp betaTemp tempFrequency idx_path