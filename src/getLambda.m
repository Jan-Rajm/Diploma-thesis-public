function LAMBDA = get_lambda(s_path, k, periodicity, nModes)
%get_lambda(s_path, nModes) calculates lambda for given k and
%   number of modes. k_x and k_y are derived from k, using variables in
%   s_patht.
%   periodicity shall be always in format p_vec = [px, py, pz], p is
%   periodicity of the examined (super)cell
dimensions = size(periodicity, 2);
switch dimensions
    case 2
        px = periodicity(1, 1);
        py = periodicity(1, 2);

        kx = s_path.xScale * k + s_path.xOffset;
        ky = s_path.yScale * k + s_path.yOffset;

        LAMBDA = diag([exp(-1j*kx*px)*ones(1, nModes),...
                       exp(-1j*ky*py)*ones(1, nModes),...
                       exp(-1j*kx*px)*ones(1, nModes),...
                       exp(-1j*ky*py)*ones(1, nModes)]);
    case 3
        px = periodicity(1, 1);
        py = periodicity(1, 2);
        pz = periodicity(1, 3);
        kx = s_path.xScale * k + s_path.xOffset;
        ky = s_path.yScale * k + s_path.yOffset;
        kz = s_path.zScale * k + s_path.zOffset;

        LAMBDA = diag([exp(-1j*kx*px)*ones(1, nModes),...
                       exp(-1j*ky*py)*ones(1, nModes),...
                       exp(-1j*kz*pz)*ones(1, nModes),...
                       exp(-1j*kx*px)*ones(1, nModes),...
                       exp(-1j*ky*py)*ones(1, nModes),...
                       exp(-1j*kz*pz)*ones(1, nModes)]);
    otherwise
        error("Unsupported number of dimensions")
end