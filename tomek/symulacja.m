clear all;

run("stale.m");

simulation_time = 2500; % Czas symulacji (sekundy)
time_step = 0.1; % Krok symulacji (sekundy)

tau_steps = tau / time_step;
tau_C_steps = tau_C / time_step;
k_min = max(tau_steps, tau_C_steps) + 1;
k_max = simulation_time / time_step + k_min;

% Trajektorie wejść procesu
T_H(1:k_max) = T_Hpp;
T_C(1:k_max) = T_Cpp;
T_D(1:k_max) = T_Dpp;
F_H(1:k_max) = F_Hpp;
F_Cin(1:k_max) = F_Cpp;
F_D(1:k_max) = F_Dpp;

F_H(500/time_step:k_max) = F_Hpp * 1.1;
F_H(1000/time_step:k_max) = F_Hpp * 0.9;
F_Cin(1500/time_step:k_max) = F_Cpp * 1.1;
F_Cin(2000/time_step:k_max) = F_Cpp * 0.9;

% Stan i wyjścia procesu przed rozpoczęciem symulacji
F(1:k_min) = alpha * sqrt(h_pp);
F_C(1:k_min) = F_Cpp;
h(1:k_min) = h_pp;
V(1:k_min) = C * h_pp^3;
T(1:k_min) = T_pp;
VT(1:k_min) = V(1) * T(1);
T_out(1:k_min) = T_pp;

% Symulacja (na podstawie modelu nieliniowego - równań różniczkowych)
for k = k_min:k_max
    dVdt = F_H(k-1) + F_C(k-1) + F_D(k-1) - F(k-1);
    dVTdt = F_H(k-1) * T_H(k-1) + F_C(k-1) * T_C(k-1) + F_D(k-1) * T_D(k-1) - F(k-1) * T(k-1);

    V(k) = V(k-1) + time_step * dVdt;
    VT(k) = VT(k-1) + time_step * dVTdt;
    T(k) = VT(k) / V(k);

    F(k) = alpha * sqrt(h(k-1));
    h(k) = nthroot(V(k) / C, 3); 
    T_out(k) = T(k - tau_steps);
    F_C(k) = F_Cin(k - tau_C_steps);
end

% figure(1);
% time = (k_min-1:k_max-1) * time_step;
% plot(time, T_out(k_min:k_max));
% hold on
% plot(time, h(k_min:k_max));
% hold off
% legend(["T_{out}", "h"]);