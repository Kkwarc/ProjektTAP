clear all;

run("stale.m");

simulation_time = 2000; % Czas symulacji (sekundy)
time_step = 1; % Krok symulacji (sekundy)

tau_steps = tau / time_step;
tau_C_steps = tau_C / time_step;
k_min = max(tau_steps, tau_C_steps) + 1;
k_max = simulation_time / time_step + k_min;

T_C(1:k_max) = T_Cpp;
T_H(1:k_max) = T_Hpp;
T_D(1:k_max) = T_Dpp;
F_Cin(1:k_max) = F_Cpp;
F_H(1:k_max) = F_Hpp;
F_D(1:k_max) = F_Dpp;

F_C(1:k_min) = F_Cpp;
F(1:k_min) = alpha * sqrt(h_pp);
h(1:k_min) = h_pp;
V(1:k_min) = C * h_pp^3;
T(1:k_min) = T_pp;
VT(1:k_min) = V(1) * T(1);
T_out(1:k_min) = T_pp;

for k = k_min:k_max
    F(k) = alpha * sqrt(h(k-1));
    F_C(k) = F_Cin(k - tau_C_steps);

    dVdt = F_H(k-1) + F_C(k-1) + F_D(k-1) - F(k-1);
    dVTdt = F_H(k-1) * T_H(k-1) + F_C(k-1) * T_C(k-1) + F_D(k-1) * T_D(k-1) - F(k-1) * T(k-1);

    V(k) = V(k-1) + time_step * dVdt;
    VT(k) = VT(k-1) + time_step * dVTdt;
    T(k) = VT(k) / V(k);
    
    h(k) = nthroot(V(k) / C, 3); 
    T_out(k) = T(k - tau_steps);
end

figure(1);
plot(T(k_min:k_max));
hold on
plot(h(k_min:k_max));
hold off
legend(["T", "h"]);



