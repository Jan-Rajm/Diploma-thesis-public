function [eig_x, eig_y] = get_eigenmode(eigenmode_file)
    eig_data = readmatrix(eigenmode_file);
    eig_data = eig_data(:, 1:2);
    eig_x = [];
    eig_y = [];


    [x, ~] = find(isnan(eig_data));
    x = unique(x);
    idxs = [];
    for i = 1:3:size(x, 1)-2
        idxs = [idxs; x(i, 1) x(i+2, 1)];
    end
    idxs = [idxs; size(eig_data, 1)+1 NaN];

    last_idx = 1;
    for i = 1:size(idxs, 1)
        curr_idx = idxs(i, 1)-1;
        size(eig_data(last_idx:curr_idx, 1));
        eig_x = [eig_x; eig_data(last_idx:curr_idx, 1)'];
        eig_y = [eig_y; eig_data(last_idx:curr_idx, 2)'];
        last_idx = idxs(i, 2)+1;
    end
end
