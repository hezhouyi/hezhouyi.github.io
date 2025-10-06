clear;
%% Parameters
% Geometry (um)
R1 = 0.5; R2 = 4;

% Diffusion coefficients (um^2/s)
D1 = 1; D2 = 2;

% Reaction, synthesis, degradation
k = 0.2;
%change k from 0.001 to 1, plot in the end 

tau = 30;  % min
mu = 1 / (tau * 60)

% Length scales
a0 = sqrt(k / D1);
a = sqrt(mu / D2);

% Partition factor
P = 1300;

% Volume and synthesis rate
V = 4* pi/3 * (R2^3 - R1^3); % so that for a cell of 4 micon, there are 100 mRNAs
N_Polysome = 0.3738*V; % so that for a cell of 4 micon, there are 100 mRNAs
lamb=0.15;
Lambda0 = lamb * N_Polysome / V;
Lambda = lamb * N_Polysome / (4*pi);

% r vectors for plotting
r1 = linspace(0.01, R1, 500);  % avoid r = 0
r2 = linspace(R1, R2, 1000);
r_all = [r1, r2];

%% Helper functions
sinh_r = @(a, r) sinh(a .* r) ./ r;
cosh_r = @(a, r) cosh(a .* r) ./ r;

% Phase I and II concentrations for Protein A
C1a_fn = @(r, C0) C0 .* sinh_r(a0, r);
C2a_fn = @(r, C1, C2, shift) C1 .* sinh_r(a, r) + C2 .* cosh_r(a, r) + shift;

% Phase I and II concentrations for Protein B
C1b_fn = @(r, A, C0) A - C0 .* sinh_r(a0, r);
C2b_fn = @(r, A1, A2) A1 .* sinh_r(a, r) + A2 .* cosh_r(a, r);

%% Common coefficients
m = (a*R2*sinh(a*R2) - cosh(a*R2)) / (a*R2*cosh(a*R2) - sinh(a*R2));
fm = a*R1 * (1 / (tanh(a*R1) - 1/m) + 1 / (1/tanh(a*R1) - m)) - 1;

%% No CTC
C0 = Lambda0*R1 / (mu * (sin(a0*R1)/P - D1 * (a0*R1*cosh(a0*R1) - sinh(a0*R1)) / (fm*D2)));
C2 = (C0*sinh(a0*R1)/P - Lambda0*R1/mu) / (cosh(a*R1) - m*sinh(a*R1));
C1 = -m * C2;

c1a_nCTC = C1a_fn(r1, C0);
c2a_nCTC = C2a_fn(r2, C1, C2, Lambda0 / mu);
c_a_nCTC = [c1a_nCTC, c2a_nCTC];

A = C0*sinh(a0*R1)/R1 + P*D1*C0*(sinh(a0*R1) - a0*R1*cosh(a0*R1))/(R1 *fm*D2);
A2 = (A*R1 - C0*sinh(a0*R1)) / (P * (cosh(a*R1) - m*sinh(a*R1)));
A1 = -m*A2;

c1b_nCTC = C1b_fn(r1, A, C0);
c2b_nCTC = C2b_fn(r2, A1, A2);
c_b_nCTC = [c1b_nCTC, c2b_nCTC];

figure(1);
A=semilogy(r1, c1a_nCTC, '-', 'LineWidth', 3,'Color',[0.4980, 0.7804, 0.9676]); hold on;
B=semilogy(r1, c1b_nCTC, 'g-', 'LineWidth', 3);

semilogy(r2, c2a_nCTC, '-', 'LineWidth', 3,'Color',[0.4980, 0.7804, 0.9676]);
semilogy(r2, c2b_nCTC, 'g-', 'LineWidth', 3); ;

xline(R1, 'k--', 'LineWidth', 1.5);
hold off;
xlabel('Radial distance from condensate center \it r \rm (\mum)');
ylabel('Protein concentration \it c \rm (\muM)');
title('No CTC Concentration Profiles');
legend([A, B], {'Protein A', 'Protein B'}, 'Location', 'northeast');
set(gca, 'FontSize', 18);
saveas(gcf, 'noCTC.svg');
writematrix([r_all.', c_a_nCTC.', c_b_nCTC.'], 'nCTC_protein_concentration_profile.csv');

% Total amount of Protein A outside the droplet (No CTC)
total_A_out_noCTC = trapz(r2, c2a_nCTC .* 4 * pi .* r2.^2);
disp(['Total material A outside (No CTC): ', num2str(total_A_out_noCTC)]);

%% CTC
n = P*(sinh(a*R1)-cosh(a*R1)/m)/sinh(a0*R1);
C1d = D2*(-sinh(a*R1)+a*R1*cosh(a*R1)+cosh(a*R1)/m-a*R1*sinh(a*R1)/m)+n*D1*(sinh(a0*R1)-a0*R1*cosh(a0*R1));
C1 = -Lambda / C1d;
C0 = n * C1;
C2 = -C1 / m;

c1a_CTC = C1a_fn(r1, C0);
c2a_CTC = C2a_fn(r2, C1, C2, 0);

A = C0*sinh(a0*R1)/R1 + P*D1*C0*(sinh(a0*R1) - a0*R1*cosh(a0*R1))/(R1 *fm*D2);
A2 = (A*R1 - C0*sinh(a0*R1)) / (P * (cosh(a*R1) - m*sinh(a*R1)));
A1 = -m*A2;

c1b_CTC = C1b_fn(r1, A, C0);
c2b_CTC = C2b_fn(r2, A1, A2);

figure(2);

A=semilogy(r1, c1a_CTC, '-', 'LineWidth', 3,'Color',[0.4980, 0.7804, 0.9676]); hold on;
B=semilogy(r1, c1b_CTC, 'g-', 'LineWidth', 3);

semilogy(r2, c2a_CTC, '-', 'LineWidth', 3,'Color',[0.4980, 0.7804, 0.9676]);
semilogy(r2, c2b_CTC, 'g-', 'LineWidth', 3); ;

xline(R1, 'k--', 'LineWidth', 1.5);
hold off;
xlabel('Radial distance from condensate center \it r \rm (\mum)');
ylabel('Protein concentration \it c \rm (\muM)');
title('CTC Concentration Profiles');
legend([A, B], {'Protein A', 'Protein B'}, 'Location', 'northeast');
set(gca, 'FontSize', 18);
saveas(gcf, 'CTC.svg');
writematrix([r_all.', c_a_nCTC.', c_b_nCTC.'], 'CTC_protein_concentration_profile.csv');

% Total amount of Protein A outside the droplet (CTC)
total_A_out_CTC = trapz(r2, c2a_CTC .* 4 * pi .* r2.^2);
disp(['Total material A outside (CTC): ', num2str(total_A_out_CTC)]);


disp(['Ratio (CTC/NoCTC): ', num2str(total_A_out_CTC / total_A_out_noCTC)]);


c_a_CTC = [c1a_CTC, c2a_CTC];  % similar to how you define c_a

% Create figure
figure(3); clf;

% First subplot: high range 
ax1 = subplot(2,1,1);
plot(r_all, c_a_nCTC, 'k-', 'LineWidth', 3); hold on;
plot(r1, c1a_CTC, '-', 'LineWidth', 3,'Color',[0.4980, 0.7804, 0.9676]); 
plot([0.5,0.5], [0,max(c1a_CTC)], '--', 'LineWidth', 3,'Color',[0.4980, 0.7804, 0.9676]);
hold off;

ylim([0.96*max(c_a_nCTC), 1.01*max(c_a_nCTC)]);
xlim([0, 4]);
set(gca, 'XTickLabel', []);
%ylabel('\it c_A \rm (\muM)');
%title('Protein A concentration profile');
legend({'No CTC','Full CTC'}, 'Location','northeast','FontSize',20);
set(gca, 'FontSize', 18);
ax1.YAxis.Exponent = 0;                    


% Second subplot: low range 
ax2 = subplot(2,1,2);
plot(r_all, c_a_nCTC, 'k-', 'LineWidth', 3); hold on;
% plot(r2, c2a_CTC, '-', 'LineWidth', 3,'Color',[0.4980, 0.7804, 0.9676]);

% Example data
x = r2;
y1 = c2a_nCTC;

% Fill the area under the curve
fill([x, fliplr(x)], [zeros(size(y1)), fliplr(y1)], ...
     [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'none');  % gray shadow
hold on;
plot(x, y1, 'k-', 'LineWidth', 3);  % Actual curve on top

y2=c2a_CTC;
fill([x, fliplr(x)], [zeros(size(y2)), fliplr(y2)], ...
     [0.3 0.3 0.3], 'FaceAlpha', 0.8, 'EdgeColor', 'none');  % gray shadow

plot([0.5,0.5], [0,max(c1a_CTC)], '--', 'LineWidth', 3,'Color',[0.4980, 0.7804, 0.9676]);
plot(r2, c2a_CTC, '-', 'LineWidth', 3,'Color',[0.4980, 0.7804, 0.9676]); hold off;
ylim([0, 1.2*max(c2a_nCTC)]);
xticks(0:0.5:10)
xlim([0, 4]);
xlabel('Radial distance from condensate center \it r \rm (\mum)');
yl = ylabel('Unmodified protein A conc. \it c_A \rm (\muM)');  % Get handle

% Get current position: [x y z]
pos = get(yl, 'Position');

% Modify position: x = a bit more left (e.g., -0.5), y = center (e.g., mean of ylim)
yl.Position = [-0.3, 1.2*max(c2a_nCTC), 0.1];

set(gca, 'FontSize', 18);

% Link x-axes for zooming/panning
linkaxes([ax1, ax2], 'x');

% Remove bottom edge of top plot
set(ax1, 'box', 'off')              % removes all box edges
set(ax1, 'XColor', 'none')          % hides x-axis (bottom) line

% Remove top edge of bottom plot
set(ax2, 'box', 'off')
set(ax2, 'XAxisLocation', 'bottom') % ensure x-axis stays on bottom
set(ax2, 'TickDir', 'out')

% Adjust spacing to make it look continuous
set(gca,'Position',[.15 .16 .8 .38])
set(ax1,'Position',[.15 .6 .8 .38])

% Save figure
saveas(gcf, 'Fig5b_A_conc_brokenY.svg');


% Define k values
k_values = logspace(-4, 2, 50);  % 50 points from 0.001 to 1
total_A_out_noCTC = zeros(size(k_values));
total_A_out_CTC = zeros(size(k_values));

for idx = 1:length(k_values)
    k = k_values(idx);

    % Update length scales
    a0 = sqrt(k / D1);

    % No CTC case -------------------------------------------
    % Common coefficients
    m = (a*R2*sinh(a*R2) - cosh(a*R2)) / (a*R2*cosh(a*R2) - sinh(a*R2));
    fm = a*R1 * (1 / (tanh(a*R1) - 1/m) + 1 / (1/tanh(a*R1) - m)) - 1;

    % Solve for concentrations
    C0 = Lambda0*R1 / (mu * (sin(a0*R1)/P - D1 * (a0*R1*cosh(a0*R1) - sinh(a0*R1)) / (fm*D2)));
    C2 = (C0*sinh(a0*R1)/P - Lambda0*R1/mu) / (cosh(a*R1) - m*sinh(a*R1));
    C1 = -m * C2;

    c2a_noCTC = C2a_fn(r2, C1, C2, Lambda0 / mu);

    total_A_out_noCTC(idx) = trapz(r2, c2a_noCTC .* 4 * pi .* r2.^2)/1e24*6.02*10^23;

    % CTC case -------------------------------------------
    n = P*(sinh(a*R1)-cosh(a*R1)/m)/sinh(a0*R1);
    C1d = D2*(-sinh(a*R1)+a*R1*cosh(a*R1)+cosh(a*R1)/m-a*R1*sinh(a*R1)/m)+n*D1*(sinh(a0*R1)-a0*R1*cosh(a0*R1));
    C1_ctc = -Lambda / C1d;
    C0_ctc = n * C1_ctc;
    C2_ctc = -C1_ctc / m;

    c2a_CTC = C2a_fn(r2, C1_ctc, C2_ctc, 0);

    total_A_out_CTC(idx) = trapz(r2, c2a_CTC .* 4 * pi .* r2.^2)/1e24*6.02*10^23;
end

% Plotting
figure(4);
loglog(k_values, total_A_out_noCTC/V, 'k-', 'LineWidth', 3);hold on;
loglog(k_values, total_A_out_CTC/V, '-', 'LineWidth', 3,'Color',[0.4980, 0.7804, 0.9676]); 
xline(0.2, 'k--', 'LineWidth', 1.5);
hold off;
xlim([10^-4,10^2])
legend({'No CTC','Full CTC'}, 'Location','northeast','FontSize',20);
xlabel('Reaction rate \it k \rm (s^{-1}) ');
ylabel({'{Mean conc. of A in dilute phase';'       (\muM)}'})
%title('Total material of A outside vs Reaction rate k');
set(gca, 'FontSize', 18,'Box','off');
saveas(gcf, 'Fig5c_cA_vs_k.svg');