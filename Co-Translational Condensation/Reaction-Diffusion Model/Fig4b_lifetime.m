format compact;
clear;

%% Parameters
R1 = 0.5;                    % Microns (fixed droplet radius)

tau1 = 0.01; tau2 = 120;      % Min lifetime range
nt = 5000;                    % Time resolution
tau = exp(linspace(log(tau1), log(tau2), nt)); % Log-spaced tau
k2_list = 1 ./ (tau * 60);   % Degradation rates in 1/s
D = 0.05;                    % Diffusion coefficient (µm^2/s)
a_list = sqrt(k2_list ./ D);
N_pro = 4 * pi * (R1 / 0.003)^3 * 0.0273 / 3;
sigma = 1.32;
NS = N_pro^(2/3) / sigma;
eps = -9.04;
v = 208 * 27 / 1e9;

%% Store 50% co-localization lifetime
tau_half_peak = [];

%% Figure 1: Selected cell sizes
selected_R2 = [2,3,4,5,6];  % Cell sizes to display curves


colors = [
    0.0500, 0.2500, 0.6000;  % Deep blue
    0.1500, 0.4000, 0.7500;  % Blue
    0.3000, 0.6000, 0.8500;  % Sky blue
    0.4980, 0.7804, 0.9676;  % tab:cyan (your reference color)
    %0.6000, 0.8500, 0.9800;  % Light cyan
    0.7000, 0.9000, 0.9900;  % Pale cyan
    %0.8000, 0.9500, 1.0000;  % Very pale cyan / close to white-blue
];

figure(1); clf; hold on;

for j = 1:length(selected_R2)
    R2 = selected_R2(j);
    V = 4/3 * pi * (R2^3 - R1^3);
    N_Polysome = 0.3738 * V;       % Number of polysomes
    frac = zeros(size(tau));
    frac_max = NS/(NS + exp(eps - 1)*V/v);
    for idx = 1:length(tau)
        a = a_list(idx);
        k2 = k2_list(idx);
        k1 = N_Polysome * k2 / V;

        m = exp(2*a*R2) * (a*R2 - 1) / (a*R2 + 1);
        n = 4 * pi * D / k2 * ((a*R1 - 1)*exp(a*R1) - m*(a*R1 + 1)*exp(-a*R1));
        l = (exp(a*R1) + m * exp(-a*R1)) / R1;
        b = k2 * l * NS - k1 * n - k2 * n * exp(eps - 1) / v;

        NB = (b - sqrt(b^2 + 4 * k1 * k2 * l * n * NS)) / (2 * k2 * l);
        frac(idx) = NB / N_Polysome;
    end
    % plot(tau, frac, 'LineWidth', 2, 'Color', colors(mod(j+3,5)+1,:), ...
    %     'DisplayName', sprintf('R_2 = %.0f \\mum', R2));
     plot(tau, frac, 'LineWidth', 2, 'Color', colors(j,:), ...
        'DisplayName', sprintf('R_2 = %.0f \\mum', R2));
    %yline(frac_max, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
end
%set(gca,"xscale","log")
% Plot decoration
yline(0.5, 'r-', 'LineWidth', 1, 'HandleVisibility', 'off');
xline(30, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

text(50, 0.54, '50% CTC', 'FontSize', 20, 'Color', 'r');
xlabel('Polysome lifetime \tau (min)');
ylabel('Fraction of CTC');
ylim([0, 1]); xlim([0, 120]);
set(gca, "FontSize", 16);
%lgd = legend('Location', 'east');
% lgd.Title.String = 'Cell size';
saveas(gcf, 'Fig4b.svg');