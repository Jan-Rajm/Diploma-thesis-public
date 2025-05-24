if simulationProfile ~= "2D_square"
   error("Structure not supported for linearized variant") 
end

m_permutation = zeros(4, 4);
m_permutation(1, 1) = 1;
m_permutation(4, 4) = 1;
m_permutation(3, 2) = 1;
m_permutation(2, 3) = 1;

m_txx = [];
m_txy = [];
m_tyx = [];
m_tyy = [];
data = sparameters("2d_conv_cst_2.s4p");%("2d_conv_cst_1.s4p");%("2d_conv_cst_1.s4p");%("s_parametry.s2p");
m_tPermutated = []; % T~

m_sParameters = data.Parameters;
frequencies = data.Frequencies;
nFrequencyPoints = size(m_sParameters, 3);
m_abcd = s2abcd(m_sParameters, data.Impedance);

px = 4.2e-3;%5.6e-3;
py = 4.2e-3;
ky = 0; % směr gama-x
kGx = ones(1, nFrequencyPoints);
kMg = ones(1, nFrequencyPoints);
k_0 = ones(1, nFrequencyPoints);
z = [];

for idx_freq = 1:nFrequencyPoints
    m_tPermutated = m_abcd(:, :, idx_freq);
    
    % rearrangement using the permutaion matrix
    m_tPermutated = m_permutation * m_tPermutated * m_permutation';
    
    m_txx = m_tPermutated(1:2, 1:2);
    m_txy = m_tPermutated(1:2, 3:4);
    m_tyx = m_tPermutated(3:4, 1:2);
    m_tyy = m_tPermutated(3:4, 3:4);
    q_y = inv(exp(-1j*ky*py)*diag([1 1]) - m_tyy);
    
    lambdaGx = eig(m_txx + m_txy * q_y * m_tyx);
    [eigVec, eigVal] = eig(m_txx + m_txy * q_y * m_tyx);
    lambdaMg = eig(m_tPermutated);
    k_0(1, idx_freq) = 2*pi*frequencies(idx_freq)/(3e8);
    
    kGx(1, idx_freq) = log(lambdaGx(1, 1))/(-1j*px);
    kMg(1, idx_freq) = log(lambdaMg(2, 1))/(-1j*px);
    
    % computation of impedance from eigenvectors
    z = [z abs(eigVec(1, 1)/eigVec(2, 2))];
end

alfaGx = abs(imag(kGx)); %útlum, čím je určeno znaménko img části?
alfaGx = alfaGx./k_0;
betaGx = abs(real(kGx)); %fázový posun

alfaMg = abs(imag(kMg));
alfaMg = alfaMg./k_0;
betaMg = abs(real(kMg));

figure
plot(alfaGx, frequencies*1e-9)
xlim([0 4e-3])
xlabel("Attenuation $\alpha/k_{0}$", "Interpreter", "latex")
ylabel("Frequency [GHz]", "Interpreter", "latex")
grid on

file_data = open("beta_2d_conv_vic.mat");
num_b = file_data.x(1, :);
num_f = file_data.x(2, :);
figure
plot(((betaGx.*(px)/pi)), frequencies*1e-9, "b-")
xlabel("Phase shift $\beta p/\pi$", "Interpreter", "latex")
ylabel("Frequency [GHz]", "Interpreter", "latex")
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex';
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';

k0 = 2 * pi * frequencies' / 3e8;

nGx = betaGx ./ k0;

figure
plot(frequencies * 1e-9, nGx)
xlabel("Frequency [GHz]", "Interpreter", "latex")
ylabel("Efective refractive index [-]", "Interpreter", "latex")
xaxisproperties = get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex';
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
ylim([1 1.5])

z0 = 120*pi;
er = nGx .* (z0 ./ z);
ur = nGx .* (z ./ z0);
plot(frequencies * 1e-9, er, frequencies * 1e-9, ur)
xlabel("Frequency [GHz]", "Interpreter", "latex")
ylabel("Relative permitivity [-], Relative permeability [-]", "Interpreter", "latex")
legend("Relative permitivity", "Relative permeability")
xaxisproperties = get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex';
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
ylim([0 5])
