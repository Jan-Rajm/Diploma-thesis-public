%% Eigenmode ASCII plot data contain also NaN at boundaries of modes. These
%  NaN values are removed and plot data is rearranged
m_eigenmodePlotData = readmatrix(addEigenmode);
eigenvalueX = [];
eigenvalueY = [];

if size(m_eigenmodePlotData, 2) >= 3
    % When more than two columns are present in raw plot data, extra
    % columns are ignored and user receives warning
    warning("Unexpected data format, additional colums ignored.")
    fprintf("Number of colums: %d\n", size(m_eigenmodePlotData, 2))
    m_eigenmodePlotData = m_eigenmodePlotData(:, 1:2);
end


[NaNposition, ~] = find(isnan(m_eigenmodePlotData));
NaNposition = unique(NaNposition);
idx_NaN = [];
for i = 1:3:size(NaNposition, 1)-2
    idx_NaN = [idx_NaN; NaNposition(i, 1) NaNposition(i+2, 1)];
end
idx_NaN = [idx_NaN; size(m_eigenmodePlotData, 1)+1 NaN];

idx_last = 1;
for i = 1:size(idx_NaN, 1)
    idx_current = idx_NaN(i, 1)-1;
    size(m_eigenmodePlotData(idx_last:idx_current, 1));
    eigenvalueX = [eigenvalueX; m_eigenmodePlotData(idx_last:idx_current, 1)'];
    eigenvalueY = [eigenvalueY; m_eigenmodePlotData(idx_last:idx_current, 2)'];
    idx_last = idx_NaN(i, 2)+1;
end

if isempty(eigenvalueX)
    warning("Unsucessfull eigenmode data accusition")
end



% clear temp idxs eig_data
