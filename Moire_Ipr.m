%% ========== band + IPR (combined plot, corrected version) ==========
clear; clc;
addpath(fullfile(fileparts(pwd), 'functions'));
%% ---- parameters ----
N1 = 12; N2 = 13;
d1 = 1;  d2 = 4;
d  = lcm(N1,N2);
k0 = 2*pi/d;
num_k =50;
rho = 1000;
c   = 1000;
% kmin=6.01533996308375
% kmax=6.01538632158324  
% kmin=6.00533996308375
% kmax=6.1  

kmin=5
kmax=7

% kmin=7.71286
% kmax=7.71332

% kmin=8.29343
% kmax=8.29889

% kmin=4.53554
% kmax=4.53724

% kmin=4.02837
% kmax=4.05544

% kmin=0
% kmax=11.1 


kmin=kmin*k0
kmax=kmax*k0
num_n = d * 8;

n1_t = gener_n(num_n, N2, d1, d2);
n2_t = gener_n(num_n, N1, d1, d2);
n_pro = n1_t + n2_t;

dz = d / num_n;

%% ---- k range ----
k_list = linspace(kmin, kmax, num_k);

%% ======================================================
% 1) band structure (run_scan)
%% ======================================================
[real_band, ~, ~] = run_scan( ...
    k_list(1), k_list(end), num_k, ...
    n_pro, d, num_n, rho, c);

%% ======================================================
% 2) IPR calculation
%% ======================================================
IPR = zeros(2, num_k);

for ii = 1:num_k

    k = k_list(ii);

    TM_T = eye(2);
    TM_series = cell(num_n,1);

    for jj = 1:num_n
        T1 = acoustic_TM(n_pro(jj), k, dz, rho, c);
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

        % normalization
        p = p ./ sqrt(sum(abs(p).^2) * dz);

        % IPR
        IPR(band,ii) = d * sum(abs(p).^4) * dz;

    end
end

%% ======================================================
% 3) A matrix (data export)
%% ======================================================
A = zeros(2*num_k, 4);

for ii = 1:num_k
    A(2*ii-1,:) = [k_list(ii)/k0, real_band(1,ii), real_band(2,ii), IPR(1,ii)];
    A(2*ii  ,:) = [k_list(ii)/k0, real_band(1,ii), real_band(2,ii), IPR(2,ii)];
end

% 可选保存
% save('band_ipr_data.mat','A');

%% ======================================================
% 4) band + IPR combined plot
%% ======================================================
figure; hold on;

scatter(k_list/k0, real_band(1,:), 18, IPR(1,:), 'filled');
scatter(k_list/k0, real_band(2,:), 18, IPR(2,:), 'filled');

colormap(turbo);
colorbar;

xlabel('k/k_0');
ylabel('qd');
grid on;
title('Band structure colored by IPR');
set(gca,'FontSize',18);



% ===== 选择最接近 k = 6.1 k0 =====

% k_target =6.01333996308375* k0;
% k_target =6.01738632158324* k0;  
%  k_target = 7.1 * k0;
k_target = 7 * k0;
% k_target = (kmax+kmin) * k0/2;
[~, idx] = min(abs(k_list - k_target));
k = k_list(idx);

TM_T = eye(2);
TM_series = cell(num_n,1);

for jj = 1:num_n
    T1 = acoustic_TM(n_pro(jj), k, dz, rho, c);
    TM_series{jj} = T1;
    TM_T = T1 * TM_T;
end

[V,~] = eig(TM_T);

band = 1;
psi_cur = V(:,band);

p = zeros(num_n,1);
z = (0:num_n-1) * dz;

for jj = 1:num_n
    p(jj) = psi_cur(1);

    if jj < num_n
        psi_cur = TM_series{jj} * psi_cur;
    end
end

% 归一化
p = p ./ sqrt(sum(abs(p).^2)*dz);

% ===== 绘图 =====
figure;
plot(z, abs(p).^2, 'LineWidth', 1.5);
xlabel('z');
ylabel('|p(z)|^2');
title(['|p(z)|^2, k/k0 = ', num2str(k/k0)]);
grid on;


B=[z'/d,abs(p).^2]
k/k0
