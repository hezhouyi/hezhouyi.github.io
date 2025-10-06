format compact;
clear;
%% Parameters
R1 = 0.5;                    % Microns (fixed droplet radius)
N_Polysome = 20;             % Number of polysomes
tau1 = 0.01; tau2 = 120;     % Lifetime range (min)
nt = 10000;                   % Resolution in tau
%tau = exp(linspace(log(tau1), log(tau2), nt)); % Log-spaced tau
tau = linspace(tau1, tau2, nt); % Linear-spaced tau

k2_list = 1 ./ (tau * 60);   % Degradation rate (1/s)
D = 0.05;                    % Diffusion coefficient (µm^2/s)
a_list = sqrt(k2_list ./ D);

N_pro = 4 * pi * (R1 / 0.003)^3 * 0.0273 / 3;
v = 208 * 27 / 1e9;

%% Parameters for plotting
R2_min = R1; R2_max = 6.5;
nR2 = 10000;
%R2_list1 = exp(linspace(log(R2_min), log(R2_max), nR2)); % log-spaced R2
R2_list1 = linspace(R2_min, R2_max, nR2);
eps_list = [-9.04, -8.97, -9.41, -8.50, -9.52];
sigma_list = [1.32, 1.30,1.54, 1.61, 1.35]; 
seq_names = ["Homo", "A-Block", "N-term", "C-term", "Cross"];
% colors = lines(length(eps_list));
colors = [
    0.4980, 0.7804, 0.9676;  % tab:cyan
    1.0000, 0.4980, 0.0549;  % tab:orange
    0.9353, 0.1686, 0.0000;  % tab:red (adjusted for MATLAB style)
    0.5961, 0.3059, 0.1235;  % tab:brown
    0.1725, 0.6275, 0.1725;  % tab:green
    0.5804, 0.4039, 0.7412   % tab:purple
];
%% Start figure
figure(2); clf; hold on;

for i = 1:length(eps_list)
    eps = eps_list(i);
    sigma = sigma_list(i);
    NS = N_pro^(2/3) / sigma;
    y = NaN(size(R2_list1));  % Store critical tau values
    
    for j = 1:length(R2_list1)
        R2 = R2_list1(j);
        V = 4/3 * pi * (R2^3 - R1^3);
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

        % Find where fraction crosses 0.5 (critical tau)
        f = find(frac(1:end-1) < 0.5 & frac(2:end) >= 0.5, 1);

        if ~isempty(f)
            y(j) = tau(f);
            R2_max(i)=R2;
        else
            y(j) = NaN;
        end
        if min(frac)>0.5
            y(j) = 0;
        end
    end

    % Plot critical tau vs. R2
    if i ==2
        plot(2*R2_list1-1, y,'--', 'LineWidth', 3, 'DisplayName', seq_names(i), 'Color', colors(i,:));
    else
        plot(2*R2_list1-1, y, 'LineWidth', 3, 'DisplayName', seq_names(i), 'Color', colors(i,:));
    end
        xline(2*R2_max(i)-1, ':', 'LineWidth', 3, 'HandleVisibility', 'off','Color', colors(i,:));
end

%% Final figure decoration
%set(gca, 'YScale', 'log');
%set(gca, 'XScale', 'log');
xlabel('Distance between condensates L (\mum)');
ylabel('Polysome life time for 50% CTC \tau (min)');
xlim([3, 11]);
xl = xlim;
xticks(floor(xl(1)):ceil(xl(2)));
ylim([0, 120]);
set(gca, 'FontSize', 16);
%legend('Location', 'northwest');
%grid on;
%title('Critical lifetime vs. Cell size');
saveas(gcf, 'Fig4d_tau_L.svg');
% set(gca, 'YScale', 'log');
% saveas(gcf, 'critical_lifetime_vs_cellsize_log.svg');