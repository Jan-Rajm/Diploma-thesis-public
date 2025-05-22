clc
% clear all
close all

data_gx = readmatrix("xmmg.txt");
data_gx = data_gx';
data_gx = data_gx(1:2, :);
x_gx = data_gx(1, 1:101);
m1_gx = data_gx(2, 1:101);
m2_gx = data_gx(2, 206:end);

data_gx1 = readmatrix("gx.txt");
data_gx1 = data_gx1';
data_xm = readmatrix("xm.txt");
data_xm = data_xm';
data_mg = readmatrix("mg.txt");
data_mg = data_mg';
x_mg = data_mg(1, 1:101);
x_i = x_mg >= 180;
m1_mg = data_mg(2, 1:101);
m2_mg = data_mg(2, 105:end);
x_gx = data_gx1(1, 1:101);
m1_gx = data_gx1(2, 1:101);
m2_gx = data_gx1(2, 105:end);

gx_idx = zeros(1, numel(x_gx));
m2_gx(1, 75:78) = [12.62 12.84 13.08 13.28];
m2_gx(1, 20:24) = linspace(14.09, 13.11, 5);
m2_gx(1, 9) = 15.12;
m2_gx(1, 47:48) = [9.35 9.15];

temp = m2_gx;
m2_gx = [m2_gx(1, 1:end/2) m1_gx(1, end/2 : end)];
m1_gx = [m1_gx(1, 1:end/2) temp(1, end/2 : end)];
m1_gx(1, 74) = 12.51;

temp = m1_mg;
m1_mg(x_i) = m2_mg(x_i);
m2_mg(x_i) = temp(x_i);


x_xm = data_xm(1, 1:(size(data_xm, 2) - 3)/2);
m1_xm = data_xm(2, 1:(size(data_xm, 2) - 3)/2);
m2_xm = data_xm(2, ((size(data_xm, 2) - 3)/2) + 4 : end);


m = [
21.6^2  21.6  1;
110^2 110 1;
136^2 136 1;
];

y = [16.24; 16.62; 16.58];

b = inv(m) * y;
b(3, 1) = b(3, 1) + 0.14;

y2 = @(x, b) b(1, 1).*x.^2 + b(2, 1).*x + b(3, 1);
x2 = 7:39;
y22 = y2(x2, b);

m2_xm(1, 7:39) = y22;
x_gx = x_gx ./ max(x_gx);
x_xm = x_xm ./ max(x_xm);
x_xm = x_xm + 1;
x_mg = x_mg ./ max(x_mg);
x_mg = x_mg + 2;
plot(x_gx, m1_gx, "bx", x_gx, m2_gx, "bx", x_xm, m1_xm,"bx", x_xm, m2_xm, "bx", x_mg, m1_mg, "bx", x_mg, m2_mg, "bx",...
    x_gx, m1_gx, "ro", x_xm, m2_xm, "ro", x_mg, m2_mg, "ro")
ax = gca;
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex';
xticks([0 1 2 3])
xticklabels(["$\Gamma$", "M", "K", "$\Gamma$"])
grid on
ylim([0 18])
xlabel("Phase shift $\beta p/\pi$", "Interpreter", "latex")
ylabel("Frequency [GHz]", "Interpreter", "latex")
% plot(x_xm, m1_xm, x_xm, m2_xm)
% f = data(:, 1)';
% n = data(:, 2)';
% 
% i1 =    1:1002;
% i2 = 1005:2006;
% i3 = 2009:3010;
% i4 = 3013:4014;
% i5 = 4016:5017;
% 
% f = f(1, i1);
% 
% plot(f, n(1, i1), f, n(1, i2), f, n(1, i3), f, n(1, i4), f, n(1, i5))
% xlabel("Frequency [GHz]", "Interpreter", "latex")
% ylabel("Effective refractive index [-]", "Interpreter", "latex")
% l = legend("$r_{i}$ = 4 mm", "$r_i$ = 3.8 mm", "$r_i$ = 2.8 mm", "$r_i$ = 2.0 mm", "$r_i$ = 0.5 mm");
% set(l, "Interpreter", "latex")

if 0
[e_x, e_y] = get_eigenmode("hex_2_GX.txt");
p_hex = 15.07e-3;
fvz = 1/((2*pi)/(40*p_hex));
% plot(e_x(1, :), e_y(1, :), e_x(1, :), e_y(2, :))
pulka_1_idx = e_x(1, :) < 180;
pulka_2_idx = e_x(1, :) > 180;
m = [e_y(1, pulka_1_idx) e_y(2, pulka_2_idx) e_y(2, pulka_1_idx) e_y(1, pulka_2_idx)];
x_xm = 1:numel(m);
x_xm = x_xm - 1;
znasobeni = 2000;
cs = [linspace(0, 20, 20), linspace(20, 0, 20)];
css = [];
mm = [];
xx = [];
x_m = 0;
for i = 1:znasobeni
    mm = [mm m];
    css = [css cs];
    xx = [xx x_xm + x_m];
    x_m = max(xx);
end
misto_nul = 1e-12;

s_cs = abs(fft(css));
s_cs = s_cs / numel(s_cs);
nuly_idx = s_cs == 0;
s_cs(1, nuly_idx) = misto_nul;
s_cs = db(s_cs)/2;
maly_idx = s_cs <= -120;
s_cs(1, maly_idx) = -120;
% figure
% plot(xx, s_cs)

s_m = abs(fft(mm));
s_m = s_m / numel(s_m);
nuly_idx = s_m == 0;
s_m(1, nuly_idx) = misto_nul;
s_m = db(s_m);
s_m = s_m / 2;
yx = 1:numel(s_m);
obalka_idx = s_m > -90;
obalka = s_m(1, obalka_idx);
obalka = [obalka; yx(1, obalka_idx)];
s_s = sin(2*pi*(1/39)*x_xm);
s_s(1, 1) = 1e-2;
s_s(1, end) = 1e-2;
t = 8*tan(pi*(x_xm/39));
l = linspace(0, fvz, numel(s_m));
figure
s_m = s_m + 90;
s_cs = s_cs + 90;
% plot(yx, s_m, obalka(2, :), obalka(1, :), yx, s_cs, "rx")
plot(l, s_m, "b-*", l, s_cs, "r--x")%, l(1, obalka_idx), obalka(1, :), "b--*",
xlim([0 fvz/2])
ylim([-0 110])
grid on
xlabel("Length~[m]", "Interpreter", "latex")
ylabel("Frequency~[$dB_{Hz}$]", "Interpreter", "latex")
legend("First mode", "Light line")
% plot(x, m, x, s_s, x, m ./ s_s, x, t)
end

if 0
eig_data = readmatrix("hex_2_GX.txt");    
eig_x = [];
eig_y = [];
[x_xm, y] = find(isnan(eig_data));
x_xm = unique(x_xm);
idxs = [];
for i = 1:3:size(x_xm, 1)-2
    idxs = [idxs; x_xm(i, 1) x_xm(i+2, 1)];
end
idxs = [idxs; size(eig_data, 1)+1 NaN];
last_idx = 1;
for i = 1:size(idxs, 1)
    curr_idx = idxs(i, 1)-1;
    size(eig_data(last_idx:curr_idx, 1));
    eig_x = [eig_x; eig_data(last_idx:curr_idx, 1)'];%+ (data_idx * 180)];
    eig_y = [eig_y; eig_data(last_idx:curr_idx, 2)'];
    last_idx = idxs(i, 2)+1;
end
% plot(eig_x(1, 1:10), eig_y(1, 1:10), eig_x(1, 11:20), eig_y(1, 11:20))
y1 = [eig_y(1, 1:10) eig_y(2, 11:20)];
y2 = [eig_y(2, 1:10) eig_y(1, 11:20)];
x1 = linspace(0, 1, 20);
x2 = linspace(1, 2, 20);
y = [y1 y2];
x_xm = [x1 x2];
% yt = 20*triangle(x/2);
plot(x_xm, y)
ylim([0 20])
end
% x = linspace(0, 3, 101);
% f1 = 0;
% f2 = 120;
% f3 = 240;
% 
% s1 = sin(2*pi*1*x + deg2rad(f1));
% s2 = sin(2*pi*1*x + deg2rad(f2));
% s3 = sin(2*pi*1*x + deg2rad(f3));
% 
% s1 = round(s1, 2);
% s2 = round(s2, 2);
% s3 = round(s3, 2);
% 
% s = [s1 s2 s3];
% 
% y = unique(s1);
% yy = unique(s);
% 
% figure
% plot(x, s1, x, s2, x, s3)
% hold on
% a = 20;
% idxs1 = s1 == y(1, a);
% idxs2 = s2 == y(1, a);
% idxs3 = s3 == y(1, a);
% plot(x(idxs1), s1(idxs1), "ro")
% hold off
% plot(linspace(0, 3, 127), yy)
% k = {};
% b = 0:0.01:1;
% for i = 1:numel(b)
%     k{1, i, 1} = [b(1, i) b(1, i)];
%     k{1, i, 2} = [abs(sin(2*pi*2*b(1, i))), abs(cos(2*pi*2*b(1, i)))];
% end
% 
% hold on
% for i = 1:101
%     plot(k{1, i, 1}, k{1, i, 2}, "bx")
% end
% hold off

% figure
% hold on
% for i = 0:0.05:1
%     u = pi*i;
%     m = 0.5;
%     a = m*exp(1j*u);
%     a_ = m*exp(-1j*u);
%     s = a + a_;
%     plot(real(a), imag(a), "rx")
%     plot(real(a_), imag(a_), "r*");
%     plot(real(s), imag(s), "bo");
%     plot([-2 2], [0 0], "k-")
%     plot([0 0], [2 -2], "k-")
%     grid on
%     xlim([-1.1 1.1])
%     ylim([-1.1 1.1])
% end
% hold off


if 0
% emproc = eigenmode_processer("hex_2_GX.txt");
emproc = eigenmode_processer("hex_super_2x2.txt");
% emproc = emproc.exclude_mode(3);
emproc = emproc.extend_modes(30);
% m1 = emproc.modes;
emproc.plot_dispersion_diagram()

emproc = emproc.sum_modes();
emproc.plot_dispersion_diagram()
emproc = emproc.calculate_spectrum(1, 8.7e-3);
emproc.plot_spectrum()
emproc = emproc.shift_spectrum_3();
% emproc = emproc.spectrum2mode();
% emproc.plot_spectrum()
% m2 = emproc.modes;
% a = max(m1.wave_nums{1, 1});
% b = max(m2.wave_nums{2, 1});
% figure
% p = plot(m1.wave_nums{1, 1}/a, m1.freqs{1, 1}, "b-", m1.wave_nums{2, 1}/a, m1.freqs{2, 1}, "b-", m2.wave_nums{2, 1}/a, m2.freqs{1, 1}*(15.4691/78.824), "r-");
% set(p,{'LineWidth'},{0.5;0.5; 0.25})
% emproc.plot_dispersion_diagram()
% emproc.plot_spectrum()
end
if 0
f = 100; % Hz
phi1 = 0;
phi2 = pi/6;
phi3 = pi/3;
phi4 = pi/2;
T = 1/f; % s
Tvz = 1e-4; 
fvz = 1/Tvz; % 
nPeriod = 10;
nSam = T * nPeriod * fvz;
A = 2;
t = linspace(0, T * nPeriod, nSam);
s1 = A * abs(sin(2*pi*f*t + phi1));
s2 = A * abs(sin(2*pi*f*t + phi2));
s3 = A * abs(sin(2*pi*f*t + phi3));
s4 = A * abs(sin(2*pi*f*t + phi4));
s = s1 .* s4;
% f_ax = fvz / nSam * (0:nSam - 1);%
f_ax = linspace(0, fvz, nSam);

figure
% plot(t, s1, "b-", t, s2, "b-", t, s3, "b-", t, s4, "b-", t, s, "r")
plot(t, s1, "b-", t, s4, "b-")
S_ = fft(s);
S = abs(fft(s))/numel(s);
% S = S(1, 1:end/2);
S1 = abs(fft(s1))/numel(s1);
% S1 = S1(1, 1:end/2);
% figure
% plot(f_ax, (real(S_)), f_ax, (imag(S_)))
plot(f_ax, S1, f_ax, S)
legend(["Původní signál", "Složený signál"])

% % [vals, idxs] = findpeaks(S);
% % idxs = round(idxs / 2, 0);
% % idxs = idxs(1, 1:10);
% % % delta_idx = mean(diff(idxs))/2;
% % % idxs = idxs - delta_idx;
% % vals = vals(1, 1:10);
% % new_spectrum = zeros(1, numel(S1));
% % new_spectrum(idxs) = vals;
% % new_spectrum(1, 1) = S(1, 1);
% % % figure
% % % plot(linspace(0, fvz/2, nSam/2), new_spectrum,"b", linspace(0, fvz/2, nSam/2), S1, "r");
% % % legend(["Nové spektrum", "Spektrum původního"])

restaurace = real(ifft(S_));
% figure
% plot(linspace(0, fvz/2, nSam), restaurace, linspace(0, fvz/2, nSam), s2)
end
% init_script_2;
% Cell_bad = Cell(1, 2);
% p_gx = Path("GX", [4 0], [0 0], Cell_2d_hex);
% p_xm = Path("XM", [0 1], [2 0], Cell_2d_hex);
% p_mg = Path("MG", [1 1], [-1 -1], Cell_2d_hex);
% p_xx = Path("xx", [0 0], [0 0], Cell_bad);
% paths = [p_gx];
% paths = [p_gx p_xx]; % zareagovalo
% f_steps = 100;
% solver = numeric_solver(1e2, 1e2, f_steps, paths, "hex_2_0_40_3m.s12p");
% ks = solver.run_mmtmm();
% figure
% hold on
% f = linspace(0, 40, f_steps);
% for i = 1:f_steps
%     plot(abs(real(ks{1, i})), f(1, i), "rx")
% end
% hold off
% % ok
% a = solver.get_impedances();
% mod = 1;
% figure
% plot(solver.selected_frequencies, real(a(mod, :)), solver.selected_frequencies, imag(a(mod, :)))
% % ok
% a = solver.get_k0();
% figure
% plot(solver.selected_frequencies, a)
% % ok
% f = [1 3 5 9];
% figure
% plot(f, solver.get_k0(f))
% T = randn(12, 12, 1);
% solver.get_determinant([1 2 3; 1j 2j 3j], T)


%% spektrum reálného a zjednodušeného modu
% clc
% close all
if 0
m1 = readmatrix("hex_2_GX_multi.txt");
m2 = m1;
m2(:, 1) = m2(:, 1) + 360;
m2(:, 2) = flip(m2(:, 2));
m = [m1; m2];
mm = [];
for i = 1:70
    nm = m;
    nm(:, 1) = nm(:, 1) + 2*(i-1)*360;
    mm = [mm; nm];
end
l = size(mm, 1);
x_xm = linspace(0, 2*i*360, l);
y = 15.4*abs(sin(2*pi*x_xm/360*(1/4)));
plot(mm(:, 1), mm(:, 2), "b-", x_xm, y)
sp = abs(fft(mm(:, 2)'))/numel(mm(:, 2));
sx = abs(fft(y))/numel(y);
figure
hold on
plot(sp, "b--")
plot(sx, "r--")
legend("opravdový", "zjednodušený")
hold off
end
% % sp(1:10) = 0;
% s = abs(ifft(sp));
% % plot(mm(:, 1)+379, s./max(s), mm(:, 1), mm(:, 2)./max(mm(:, 2)))
%%










%%
% mmtmm_main_zaloha;
% b = ALL_B{1, 1};
% b = round(b, 4);
% f = ALL_F{1, 1};
% f = round(f, 4);
% f = f*1e-9;
% f = f/max(f);
% 
% y = 0.6352*abs(sin(b/(pi/4)));
% 
% plot(b, f, "bx", b, y, "ro")

%  GX
% f1 = ALL_F{1, 1}/max(ALL_F{1, 1});   
% f2 = ALL_F{2, 1}/max(ALL_F{2, 1});   
% x = linspace(0, 2, 40);
% y1 = 0.37*abs(cos(2*pi*x*(1/8) + (pi/2) ));
% y2 = 0.18*    cos(2*pi*x*(1/4) + 0*pi)    + 0.72;
% y3 = 0.40*   (cos(2*pi*x*(1/2)).^2)       + 0.95-0.35;

% % %  XM
% % % f1 = ALL_F{1, 1}/max(ALL_F{1, 1});   
% % % f2 = ALL_F{2, 1}/max(ALL_F{2, 1});   
% % % x = linspace(0, 2, 40);
% % % y1 = 0.56*abs(sin(2*pi*x*(1/8) + (pi/2) ));
% % % y2 = 0.21*    cos(2*pi*x*(1/4) + 1*pi)    + 0.58;
% % % y3 = 0.40*   (cos(2*pi*x*(1/2)).^2)       + 0.95-0.35;

%  plot(ALL_B{1, 1}*2, f1, "bx", ALL_B{2, 1}*2, f2, "bx", x, y1, "ro", x, y2, "go")%, x, y3, "mo")

% 
% if f(b == min(b)) == min(f) % když počáteční bod je úplně vlevo dole
%     x_ext = [min(b)];
%     y_ext = [min(f)];
%     b = b(1, 2:end);
%     f = f(1, 2:end);
% else
%     error("Nefunguje")
% end
% x_ext = [0.0404 0.0707 0.0909 0.1111];
% y_ext = [0.8098 1.2097 1.6096 2.0494];

%% nalezení dalších tří bodů pro extrapolaci
% for pts_processed = 1:3
%     min_distance = Inf;
%     n_pts = 1;
%     while(numel(b) > n_pts)
%         current_pt_x = b(1, n_pts);
%         current_pt_y = f(1, n_pts);
%         distance = get_hypotenuse(x_ext(1, end) - current_pt_x, y_ext(1, end) - current_pt_y);
%         if distance < min_distance && distance > 0
%             min_distance = distance;
%             pt_to_add = [current_pt_x current_pt_y];
%             idx_to_pop = (b == current_pt_x) & (f == current_pt_y);
%             b(idx_to_pop) = NaN;
% %             b = b(~idx_to_pop);
% %             f = f(~idx_to_pop);
%         end
%         n_pts = n_pts + 1;
%     end
%     x_ext = [x_ext pt_to_add(1, 1)];
%     y_ext = [y_ext pt_to_add(1, 2)];
% end
% % ext = interp1(x, y, delimiters(path+1, 1), "pchip", "extrap")
% x_unique = sort(unique(b(b > max(x_ext))));
% x_unique = x_unique(x_unique > min(b));
% x_unique = x_unique(~isnan(x_unique));
% %% extrapolace, hledání dalších padnoucích bodů
% for next_x = x_unique
%     new_y = interp1(x_ext, y_ext, next_x, "pchip", "extrap");
%     distance = Inf;
%     for pts_to_compare = b(b > max(x_ext))
%         
%     end
%     x_ext = [x_ext next_x];
%     y_ext = [y_ext new_y];
% end
% 
% plot(b, f, "bx", x_ext, y_ext, "ro", b, f, "g*")
% ylim([0 1])



% b = round(ALL_B{1, 1}, 3);
% b_ = b(b < 0.5);
% f_ = f(b < 0.5);
% f_ = f_(f_ < 0.8);
% b_ = b_(f_ < 0.8);
% b__ = flip_plot(b_)/max(b__)+0.5;
% 
% f = ALL_F{1, 1};
% f = f*1e-9;
% max_f = max(f);
% f = f/max_f;
% 
% 
% plot(b_, f_, b__, f_)
% 
% [B1, B2] = meshgrid(b, b);
% [F1, F2] = meshgrid(f, f);
% 
% B_MAT = abs(B1 - B2);
% F_MAT = abs(F1 - F2);

% close all
% clc
% clear b f x_ext y_ext 
% 

% y = f;
% b_1 = b(b < 0.5);
% f_1 = f(b < 0.5);
% 
% x_ext = [];
% y_ext = [];
% 
% if f(b == min(b)) == min(f) % když počáteční bod je úplně vlevo dole
%     x_ext = [min(b)];
%     y_ext = [min(f)];
%     b = b(1, 2:end);
%     f = f(1, 2:end);
% else
%     error("Nefunguje")
% end
% %% nalezení dalších tří bodů pro extrapolaci
% for pts_processed = 1:3
%     min_distance = Inf;
%     n_pts = 1;
%     while(numel(b) > n_pts)
%         current_pt_x = b(1, n_pts);
%         current_pt_y = f(1, n_pts);
%         distance = get_hypotenuse(x_ext(1, end) - current_pt_x, y_ext(1, end) - current_pt_y);
%         if distance < min_distance && distance > 0
%             min_distance = distance;
%             pt_to_add = [current_pt_x current_pt_y];
%             idx_to_pop = (b == current_pt_x) & (f == current_pt_y);
%             b(idx_to_pop) = NaN;
% %             b = b(~idx_to_pop);
% %             f = f(~idx_to_pop);
%         end
%         n_pts = n_pts + 1;
%     end
%     x_ext = [x_ext pt_to_add(1, 1)];
%     y_ext = [y_ext pt_to_add(1, 2)];
% end
% % ext = interp1(x, y, delimiters(path+1, 1), "pchip", "extrap")
% x_unique = sort(unique(b(b > max(x_ext))));
% x_unique = x_unique(x_unique > min(b));
% x_unique = x_unique(~isnan(x_unique));
% %% extrapolace, hledání dalších padnoucích bodů
% for next_x = x_unique
%     new_y = interp1(x_ext, y_ext, next_x, "pchip", "extrap");
%     distance = Inf;
%     for pts_to_compare = b(b > max(x_ext))
%         
%     end
%     x_ext = [x_ext next_x];
%     y_ext = [y_ext new_y];
% end
% 
% plot(x, y, "bx", x_ext, y_ext, "ro", b, f, "g*")

% 
% 
% 
% 
% 
% 











% b = round(ALL_B, 3);
% f = ALL_F(b ~= 0.5);
% f = f*1e-9;
% b = b(b ~= 0.5);
% 
% 
% hold on 
% plot(b, f, "x")
% % ext = interp1(x, y, delimiters(path+1, 1), "pchip", "extrap")
% f_span = unique(f);
% f_span = f_span(1, 1:4);
% b_unique = sort(unique(b));
% k = [];
% for i = 1:4
%     k = [k min(b(f == f_span(i)))];
% end
% 
% for i = 5:size(b_unique, 2)
% %     next_b
% %     next_f = interp1(k, f_span, )
% end

% hold off
% clear b f k
% close all
% clc
% 
% add_eigenmode = "hex_2_GX_XM_MG.txt";
% process_eigenmode;
% 
% eig_y1 = eig_y;
% eig_y2 = eig_y;
% 
% eig_x1 = eig_x;
% 
% eig_x2 = flip_plot(eig_x);
% eig_x2 = eig_x2 + 1;
% for i = 1:4
%     eig_x2(i, :) = flip(eig_x2(i, :));
%     eig_y2(i, :) = flip(eig_y2(i, :));
% end
% eig_a = eig_y1(1, :);
% eig_a(1, 1:end-5) = NaN;
% eig_b = eig_y2(2, :);
% eig_b(2, 5:end) = NaN;
% %% Krok 1
% hold on
% for i = 1:4
%     plot(eig_x1(i, :), eig_y1(i, :), "b-")
% end
% plot(eig_x1(1, :), eig_a, "rx")
% 
% for i = 1:4
%     plot(eig_x2(i, :), eig_y2(i, :), "b-")
% end
% plot(eig_x2(1, :), eig_b, "ro")
% hold off
% 
% temp = eig_y2(2, :);
% eig_y2(2, :) = eig_y2(1, :);
% eig_y2(1, :) = temp;
% 
% temp = eig_y2(3, :);
% eig_y2(3, :) = eig_y2(2, :);
% eig_y2(2, :) = temp;
% 
% %% Krok 2
% figure
% hold on
% for i = 1:4
%     for ii = 1:4
%         if and(sum(sign(diff([eig_y1(i, end-5:end) eig_y2(ii, 1:5)]))) > 8,...
%                 round(eig_y(i, end) - eig_y2(ii, 1)) == 0)
%             sum(diff([eig_y1(i, end-5:end) eig_y2(ii, 1:5)]))
%             plot(eig_x1(i, :), eig_y1(i, :))
%             plot(eig_x2(ii, :), eig_y2(ii, :))
%             
%         end
% %         sum(sign(diff([eig_y1(i, end-5:end) eig_y2(i, 1:5)])))
%     end
% end
% hold off