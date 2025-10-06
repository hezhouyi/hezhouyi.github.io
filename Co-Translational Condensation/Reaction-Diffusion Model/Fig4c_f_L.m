format compact;
clear;
%% Parameters
R1 = 0.5;                    % Microns (fixed droplet radius)
tau = 30; % min

k2= 1 / (tau * 60);   % Degradation rate (1/s)
D = 0.05;                    % Diffusion coefficient (µm^2/s)
a= sqrt(k2/ D);

N_pro = 4 * pi * (R1 / 0.003)^3 * 0.0273 / 3;

v = 208 * 27 / 1e9;

%% Parameters for plotting
R2_list1 = linspace(0.5,6.5, 1000);


% eps_list = [-9.04, -8.97, -9.41, -8.50, -9.53];
% sigma_list = [1.32, 1.30,1.54, 1.61, 1.35];
eps_list = [-9.07, -8.99, -9.44, -8.50, -9.57];
sigma_list = [1.33, 1.31,1.56, 1.61, 1.36];
seq_names = ["Homo", "M-Block", "N-term", "C-term", "Cross"];
colors = [
    0.4980, 0.7804, 0.9676;  % tab:cyan
    1.0000, 0.4980, 0.0549;  % tab:orange
    0.9353, 0.1686, 0.0000;  % tab:red (adjusted for MATLAB style)
    0.5961, 0.3059, 0.1235;  % tab:brown
    0.1725, 0.6275, 0.1725;  % tab:green
    0.5804, 0.4039, 0.7412   % tab:purple
];
%% Start figure
fig=figure(1); clf; hold on;
set(gcf, 'Position', [584 495 560 420]);  % [x, y, width, height]
for i = 1:length(eps_list)
    eps = eps_list(i);
    sigma = sigma_list(i);
    NS = N_pro^(2/3) / sigma;
    y = NaN(size(R2_list1));  % Store critical tau values
    
    for j = 1:length(R2_list1)
        R2 = R2_list1(j);
        V = 4/3 * pi * (R2^3 - R1^3);
            N_Polysome = 0.3738*V;             % Number of polysomes
            k1 = N_Polysome * k2 / V;
            m = exp(2*a*R2) * (a*R2 - 1) / (a*R2 + 1);
            n = 4 * pi * D / k2 * ((a*R1 - 1)*exp(a*R1) - m*(a*R1 + 1)*exp(-a*R1));
            l = (exp(a*R1) + m * exp(-a*R1)) / R1;
            b = k2 * l * NS - k1 * n - k2 * n * exp(eps - 1) / v;

            NB = (b - sqrt(b^2 + 4 * k1 * k2 * l * n * NS)) / (2 * k2 * l);
            frac(i,j) = NB / N_Polysome;

    end

    %Plot critical tau vs. R2
    if i ==2
        plot(2*R2_list1-1, frac(i,:),'--', 'LineWidth', 3, 'DisplayName', seq_names(i), 'Color', colors(i,:));
    else
        plot(2*R2_list1-1, frac(i,:), 'LineWidth', 3, 'DisplayName', seq_names(i), 'Color', colors(i,:));
    end
end

%% Final figure decoration

xlabel('Distance between condensates L (\mum)');
ylabel('Fraction of CTC    ');
xlim([0, 12]);
xl = xlim;
xticks(floor(xl(1)):ceil(xl(2)));
ylim([0, 1]);
set(gca, 'FontSize', 16);
legend('Location', 'northeast', 'FontSize', 24);
% Save figure
set(gca, 'LooseInset', get(gca, 'TightInset'));
print(fig, 'Fig4c_f_L.svg', '-dsvg');

% Save data
frac(:,1)=1;
writematrix(frac', 'frac30.csv');