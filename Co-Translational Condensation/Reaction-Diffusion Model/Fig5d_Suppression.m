format compact;
clear;

%% Parameters
R1 = 0.5;                    % Microns (fixed droplet radius)
R_cell = 10;
V_cell = 4/3 * pi * (R_cell^3 - R1^3); % constant of cell volume

R2_list1 = linspace(0.5, 6.5, 1000);

f = readmatrix('frac30.csv');  % f of 30 min life time
f = f';


% Diffusion coefficients (um^2/s)
D1 = 1; D2 = 2;

% Reaction, synthesis, degradation
k = 0.5;
tau = 30;  % min
mu = 1 / (tau * 60);

lamb =0.15;
a0 = sqrt(k / D1);
a = sqrt(mu / D2);

P = 1200;

N_pro = 4 * pi * (R1 / 0.003)^3 * 0.0273 / 3;
v = 208 * 27 / 1e9;

eps_list = [-9.04, -8.97, -9.41, -8.50, -9.53];
sigma_list = [1.32, 1.30,1.54, 1.61, 1.35];
seq_names = ["Homo", "M-Block", "N-term", "C-term", "Cross"];

colors = [
    0.4980, 0.7804, 0.9676;  % tab:cyan
    1.0000, 0.4980, 0.0549;  % tab:orange
    0.9353, 0.1686, 0.0000;  % tab:red
    0.5961, 0.3059, 0.1235;  % tab:brown
    0.1725, 0.6275, 0.1725;  % tab:green
    0.5804, 0.4039, 0.7412   % tab:purple
];

%% Helper functions
sinh_r = @(a, r) sinh(a .* r) ./ r;
cosh_r = @(a, r) cosh(a .* r) ./ r;

% Phase I and II concentrations for Protein A
C1a_fn = @(r, C0) C0 .* sinh_r(a0, r);
C2a_fn = @(r, C1, C2, shift) C1 .* sinh_r(a, r) + C2 .* cosh_r(a, r) + shift;

% Phase I and II concentrations for Protein B
C1b_fn = @(r, A, C0) A - C0 .* sinh_r(a0, r);
C2b_fn = @(r, A1, A2) A1 .* sinh_r(a, r) + A2 .* cosh_r(a, r);

%% Initialize storage arrays
N_Dilute_noCTC = zeros(length(eps_list), length(R2_list1));
N_Dilute_CTC = zeros(length(eps_list), length(R2_list1));
c_Dilute_noCTC = zeros(length(eps_list), length(R2_list1));
c_Dilute_CTC = zeros(length(eps_list), length(R2_list1));
E = zeros(length(eps_list), length(R2_list1));
S = zeros(length(eps_list), length(R2_list1));

%% Main calculation loop
for i = 1:length(eps_list)
    eps = eps_list(i);
    sigma = sigma_list(i);
    NS = N_pro^(2/3) / sigma;
    for j = 1:length(R2_list1)
        R2 = R2_list1(j);
        V = 4/3 * pi * (R2^3 - R1^3); % shell of dilute phase with one droplet
        N_Polysome = 0.3738 * V;       % Number of polysomes
        N_D = V_cell / V;
        
        Lambda0 = lamb * N_Polysome / V;
        Lambda = lamb * N_Polysome / (4*pi);
        
        r2 = linspace(R1, R2, 1000);
        
        %% Common coefficients
        m = (a*R2*sinh(a*R2) - cosh(a*R2)) / (a*R2*cosh(a*R2) - sinh(a*R2));
        fm = a*R1 * (1 / (tanh(a*R1) - 1/m) + 1 / (1/tanh(a*R1) - m)) - 1;
        
        %% No CTC calculation
        C0 = Lambda0*R1 / (mu * (sin(a0*R1)/P - D1 * (a0*R1*cosh(a0*R1) - sinh(a0*R1)) / (fm*D2)));
        C2 = (C0*sinh(a0*R1)/P - Lambda0*R1/mu) / (cosh(a*R1) - m*sinh(a*R1));
        C1 = -m * C2;
        
        c2a_nCTC = C2a_fn(r2, C1, C2, Lambda0 / mu);
        N_Dilute_noCTC(i,j) = trapz(r2, c2a_nCTC .* 4 * pi .* r2.^2);
        c_Dilute_noCTC(i,j) = N_Dilute_noCTC(i,j) / V;
        
        %% CTC calculation
        n = P*(sinh(a*R1)-cosh(a*R1)/m)/sinh(a0*R1);
        C1d = D2*(-sinh(a*R1)+a*R1*cosh(a*R1)+cosh(a*R1)/m-a*R1*sinh(a*R1)/m) + ...
              n*D1*(sinh(a0*R1)-a0*R1*cosh(a0*R1));
        C1 = -Lambda / C1d;
        C0 = n * C1;
        C2 = -C1 / m;
        
        c2a_CTC = C2a_fn(r2, C1, C2, 0);
        N_Dilute_CTC(i,j) = trapz(r2, c2a_CTC .* 4 * pi .* r2.^2);
        c_Dilute_CTC(i,j) = N_Dilute_CTC(i,j) / V;
        
        %% Calculate efficiency and surplus
        E(i,j) = f(i,j) * c_Dilute_CTC(i,j) / c_Dilute_noCTC(i,j) + (1 - f(i,j));
        S(i,j) = f(i,j) * c_Dilute_CTC(i,j) + (1 - f(i,j)) * c_Dilute_noCTC(i,j);
    end
end
E(:,1)=1;

% Combined figure
fig1 = create_combined_figure(R2_list1, S, E, c_Dilute_noCTC, seq_names, colors);
set(fig1, 'PaperPositionMode', 'auto');
print(fig1, 'Fig5d_Suppression.svg', '-dsvg', '-painters');

function fig = create_combined_figure(R2_list, S, E, c_noCTC, seq_names, colors)
    fig = figure('Name', 'Combined Analysis');
    
    % Top subplot - Surplus (log scale)
    ax1 = subplot(2, 1, 1);
    hold on;
    
    for i = 1:size(S, 1)
        if i == 2  % M-Block sequence gets dashed line
            plot(2*R2_list-1, S(i,:), '--', 'LineWidth', 2.5, ...
                 'DisplayName', seq_names(i), 'Color', colors(i,:));
        else
            plot(2*R2_list-1, S(i,:), '-', 'LineWidth', 2.5, ...
                 'DisplayName', seq_names(i), 'Color', colors(i,:));
        end
    end
    % Add no CTC reference line
    plot(2*R2_list-1, c_noCTC(end,:), 'k-', 'LineWidth', 2.5, 'DisplayName', 'No CTC');
    
    set(gca, 'YScale', 'log');
    ylim([10^-4.5, 10^0.5]);
    xlim([0, 10]);
    % Explicitly set ticks
    yticks([1e-4,1e-3, 1e-2, 1e-1, 1e0]);  % 10^-3 to 10^0
    ylabel({'Mean conc. of A','      (\muM)'});
    set(gca, 'FontSize', 16);
    set(gca, 'Position', [0.16, 0.58, 0.8, 0.38]);
    
    legend('Location', 'southeast', 'FontSize', 16);
    
    % Bottom subplot - Efficiency
    ax2 = subplot(2, 1, 2);
    hold on;
    
    for i = 1:size(E, 1)
        if i == 2  % M-Block sequence gets dashed line
            plot(2*R2_list-1, E(i,:), '--', 'LineWidth', 3, ...
                 'DisplayName', seq_names(i), 'Color', colors(i,:));
        else
            plot(2*R2_list-1, E(i,:), '-', 'LineWidth', 3, ...
                 'DisplayName', seq_names(i), 'Color', colors(i,:));
        end
    end
    
    yline(1, 'k-', 'LineWidth', 2.5);

    
    xlabel('Distance between condensates L (\mum)');
    ylabel({'Relative surplus'; '   '})
    % ylabel({'Relative surplus'; '<c_A^{i}>/<c_A^{No CTC}>'})
    ylim([0, 1]);
    xlim([0, 10]);
    set(gca, 'FontSize', 16);
    set(gca, 'Position', [0.16, 0.13, 0.8, 0.38]);


end