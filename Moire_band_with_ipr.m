% Moiré band structure + IPR (coarse + fine, NO run_scan modification)
clear; clc; close all;

script_dir = fileparts(mfilename('fullpath'));
figure_code_dir = fileparts(script_dir);
addpath(fullfile(figure_code_dir, 'functions'));

%% ========== Parameters ==========
N1 = 12;
N2 = 13;
d1 = 1;
d2 = 4;

d  = lcm(N1,N2);
k0 = 2*pi/d;

k_min0 = 0;
k_max0 = 11.1*k0;

num_k_coarse = 1001;
num_k_fine   = 200;

rho = 1000;
c   = 1000;

slope_threshold = pi;
width_threshold = 0.5;
amplitude_threshold = 1;

max_refinement_iter = 25;
refinement_margin = 0.3;
refinement_points_factor = 1.5;

%% ========== Narrow band search ==========
fprintf('Finding narrow band...\n');

[final_result, refinement_history, candidates, k_coarse, real_coarse] = ...
    find_narrow_band( ...
    N1, N2, d1, d2, ...
    k_min0, k_max0, ...
    num_k_coarse, num_k_fine, ...
    slope_threshold, width_threshold, amplitude_threshold, ...
    max_refinement_iter, refinement_margin, refinement_points_factor, ...
    rho, c, true);

%% ========== lattice ==========
num_n = d * 8;

n1_t = gener_n(num_n, N2, d1, d2);
n2_t = gener_n(num_n, N1, d1, d2);
n_pro = n1_t + n2_t;

dz = d / num_n;

%% =========================================================
%                  COARSE IPR (IMPORTANT)
%% =========================================================
fprintf('Computing coarse IPR...\n');

IPR_coarse = zeros(2, length(k_coarse));

for ii = 1:length(k_coarse)

    k = k_coarse(ii);

    TM_T = eye(2);
    TM_series = cell(num_n,1);

    for jj = 1:num_n
        S_curr = n_pro(jj);
        T1 = acoustic_TM(S_curr, k, dz, rho, c);

        TM_series{jj} = T1;
        TM_T = T1 * TM_T;
    end

    [V,~] = eig(TM_T);

    for band = 1:2

        psi = V(:,band);

        psi_cur = psi;
        p = zeros(num_n,1);

        for jj = 1:num_n
            p(jj) = psi_cur(1);
            if jj < num_n
                psi_cur = TM_series{jj} * psi_cur;
            end
        end

        p = p ./ sqrt(sum(abs(p).^2) * dz);

        IPR_coarse(band,ii) = d * sum(abs(p).^4) * dz;

    end
end

%% ========== fine scan (dispersion only from run_scan) ==========
[real_final, imag_final, k_final] = run_scan( ...
    final_result(1), final_result(2), num_k_fine*2, ...
    n_pro, d, num_n, rho, c);

%% =========================================================
%                  FINE IPR (IMPORTANT)
%% =========================================================
fprintf('Computing fine IPR...\n');

IPR_final = zeros(2, length(k_final));

for ii = 1:length(k_final)

    k = k_final(ii);

    TM_T = eye(2);
    TM_series = cell(num_n,1);

    for jj = 1:num_n
        S_curr = n_pro(jj);
        T1 = acoustic_TM(S_curr, k, dz, rho, c);

        TM_series{jj} = T1;
        TM_T = T1 * TM_T;
    end

    [V,~] = eig(TM_T);

    for band = 1:2

        psi = V(:,band);

        psi_cur = psi;
        p = zeros(num_n,1);

        for jj = 1:num_n
            p(jj) = psi_cur(1);
            if jj < num_n
                psi_cur = TM_series{jj} * psi_cur;
            end
        end

        p = p ./ sqrt(sum(abs(p).^2) * dz);

        IPR_final(band,ii) = d * sum(abs(p).^4) * dz;

    end
end

%% ===================== PLOT 1: COARSE =====================
figure;
hold on;

scatter(k_coarse/k0, real_coarse(1,:), 12, IPR_coarse(1,:), 'filled');
scatter(k_coarse/k0, real_coarse(2,:), 12, IPR_coarse(2,:), 'filled');

colormap(turbo);
colorbar;

xlabel('k/k_0');
ylabel('qd');
title('Coarse band structure with IPR');
grid on;

%% ===================== PLOT 2: FINE =====================
figure;
hold on;

scatter(k_final/k0, real_final(1,:), 35, IPR_final(1,:), 'filled');
scatter(k_final/k0, real_final(2,:), 35, IPR_final(2,:), 'filled');

colormap(turbo);
colorbar;

xlabel('k/k_0');
ylabel('qd');
title('Fine band structure with IPR');
grid on;