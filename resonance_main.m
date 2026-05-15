%% Lab 3: Resonance - Main Driver Script
% Calls idealcase_sim.m and nonidealcase_sim.m


clear; clc; close all;


%  SHARED SETTINGS

v0  = 0;       % Initial condition: v(0)
u0  = 0;       % Initial condition: u(0) = v'(0)
tol = 1e-6;    % Convergence tolerance for periodic orbit

%  PART I: Ideal Case - Amplitude Response Curve
%  gamma = 0.02, omega swept from 0.5 to 1.5

fprintf('=== PART I: Ideal Case ===\n');

gamma_I    = 0.02;
omega_vec  = linspace(0.5, 1.5, 501);
R_numeric  = zeros(size(omega_vec));
R_analytic = 1 ./ sqrt((1 - omega_vec.^2).^2 + 4*gamma_I^2*omega_vec.^2);

% Warm-start: carry converged state forward to next omega
v0_warm = v0;
u0_warm = u0;

for k = 1:length(omega_vec)
    omega = omega_vec(k);
    [~, v_tmp, u_tmp, R_numeric(k)] = idealcase_sim(v0_warm, u0_warm, gamma_I, omega, tol);
    v0_warm = v_tmp(end);
    u0_warm = u_tmp(end);

    if mod(k-1, 10) == 0
        fprintf('  omega = %.2f | R_numeric = %.4f | R_analytic = %.4f\n', ...
                omega, R_numeric(k), R_analytic(k));
    end
end

figure('Name', 'Part I: Ideal Amplitude Response');
plot(omega_vec, R_analytic, 'k-',  'LineWidth', 2,   'DisplayName', 'Analytic'); hold on;
plot(omega_vec, R_numeric,  'r--', 'LineWidth', 1.5, 'DisplayName', 'Numeric');
xlabel('\omega'); ylabel('Response Amplitude R');
title('Part I: Ideal Case Amplitude Response (\gamma = 0.02)');
legend; grid on;

%  PART II: Nonideal Case - Bistability
%  epsilon = 0.001, gamma = 0.06 / 0.04 / 0.02
%  Sweep omega both UP (0.5->1.5) and DOWN (1.5->0.5)

fprintf('\n=== PART II: Nonideal Case ===\n');

epsilon = 0.001;
gamma_vals = [0.06, 0.04, 0.02];
omega_up = linspace(0.5, 1.5, 500);
omega_dn = linspace(1.5, 0.5, 500);
colors = lines(length(gamma_vals)*2 + 1);

for gi = 1:length(gamma_vals)
    gamma = gamma_vals(gi);
    R_up = zeros(size(omega_up));
    R_dn = zeros(size(omega_dn));

    % Sweep up in omega
    v0_warm = v0; u0_warm = u0;
    for k = 1:length(omega_up)
        omega = omega_up(k);
        [~, v_tmp, u_tmp, R_up(k)] = nonidealcase_sim(v0_warm, u0_warm, gamma, omega, epsilon, tol);
        v0_warm = v_tmp(end);
        u0_warm = u_tmp(end);
    end

    % Sweep down in omega
    v0_warm = v0; u0_warm = u0;
    for k = 1:length(omega_dn)
        omega = omega_dn(k);
        [~, v_tmp, u_tmp, R_dn(k)] = nonidealcase_sim(v0_warm, u0_warm, gamma, omega, epsilon, tol);
        v0_warm = v_tmp(end);
        u0_warm = u_tmp(end);
    end

    fprintf('  gamma = %.2f done\n', gamma);
    fprintf('      R up max = %.2f, occurs at omega = %.2f \n', max(R_up), omega_up(R_up == max(R_up)));
    fprintf('      R dn max = %.2f, occurs at omega = %.2f \n', max(R_dn), omega_dn(R_dn == max(R_dn)));

    figure('Name', 'Part II: Nonideal Amplitude Response - Bistability');

    plot(omega_up, R_up, 'Color', colors(gi*2,:), 'LineWidth', 2, ...
         'DisplayName', 'omega sweeping up'); hold on;
    plot(omega_dn, R_dn, 'Color', colors(gi*2 + 1,:), 'LineWidth', 1.5, ...
         'DisplayName', 'omega sweeping down');

    xlabel('\omega'); ylabel('Response Amplitude R');
    title(sprintf('Nonideal Case - Bistability ($\\epsilon$ = 0.001)($\\gamma$ = %.3f)', gamma), ...
        'Interpreter', 'latex');
    legend('Location', 'northwest'); grid on;

end

%  PART III: Strongly Nonlinear - Phase Portrait + Poincare Section
%  gamma = 0.09, omega = 0.8, epsilon = 800

fprintf('\n=== PART III: Strongly Nonlinear Case ===\n');

gamma_III = 0.09;
omega_III = 0.8;
epsilon_III = 800;
n_periods = 5000;

% Converge to periodic orbit first
[~, v_conv, u_conv, ~] = nonidealcase_sim(v0, u0, gamma_III, omega_III, epsilon_III, tol);

% Phase portrait of converged orbit
figure('Name', 'Part III: Phase Portrait');
plot(v_conv, u_conv, 'b-', 'LineWidth', 0.3);
xlabel('v'); ylabel('dv/dt');
title('Part III: Phase Portrait (\gamma=0.09, \omega=0.8, \epsilon=800, periods = 200)');
grid on;

% Build Poincare section: sample [v, u] at start of each period
poincare_v = zeros(n_periods, 1);
poincare_u = zeros(n_periods, 1);
v0_p = v_conv(end);
u0_p = u_conv(end);

for k = 1:n_periods
    [~, v_tmp, u_tmp, ~] = nonidealcase_sim(v0_p, u0_p, gamma_III, omega_III, epsilon_III, NaN);
    poincare_v(k) = v_tmp(end);
    poincare_u(k) = u_tmp(end);
    v0_p = v_tmp(end);
    u0_p = u_tmp(end);
end

figure('Name', 'Part III: Poincare Section');
plot(poincare_v, poincare_u, 'k.', 'MarkerSize', 4);
xlabel('v'); ylabel('dv/dt');
title(sprintf(['Part III: Poincar\\''{e} Section (%d periods, ', ...
               '$\\gamma=0.09$, $\\omega=0.8$, $\\varepsilon=800$)'], ...
               n_periods), ...
      'Interpreter', 'latex');
grid on;

%  BONUS: Find another interesting case
%  gamma=0.05, omega=0.85, epsilon=800 (period-doubling candidate)

fprintf('\n=== BONUS: Alternative Parameter Exploration ===\n');

% parameters (chosen by hand after some running of the code to find which produce an interesting result)
gamma_b   = 1;
omega_b   = 0.02;
epsilon_b = 1100;

% Converge to periodic oscillation
[~, v_b, u_b, ~] = nonidealcase_sim(v0, u0, gamma_b, omega_b, epsilon_b, tol);

% Phase portrait of converged orbit
figure('Name', 'Bonus: Phase Portrait');
plot(v_b, u_b, 'b-', 'LineWidth', 0.3);
xlabel('v'); ylabel('dv/dt');
title('Part III: Phase Portrait (\gamma=1, \omega=0.02, \epsilon=1100, periods = 200)');
grid on;

% Build Poincare section
poincare_vb = zeros(n_periods, 1);
poincare_ub = zeros(n_periods, 1);
v0_b = v_b(end);
u0_b = u_b(end);

for k = 1:n_periods
    [~, v_tmp, u_tmp, ~] = nonidealcase_sim(v0_b, u0_b, gamma_b, omega_b, epsilon_b, NaN);
    poincare_vb(k) = v_tmp(end);
    poincare_ub(k) = u_tmp(end);
    v0_b = v_tmp(end);
    u0_b = u_tmp(end);
end

figure('Name', 'Bonus: Poincare Section');
plot(poincare_vb, poincare_ub, 'r.', 'MarkerSize', 4);
xlabel('v'); ylabel('dv/dt');
title(sprintf("Bonus Poincar\\'{e} Section ($\\gamma$=%.2f, $\\omega$=%.2f, $\\epsilon$=%d)", ...
              gamma_b, omega_b, epsilon_b), ...
      'Interpreter', 'latex');
grid on;

fprintf('\nAll parts complete.\n');