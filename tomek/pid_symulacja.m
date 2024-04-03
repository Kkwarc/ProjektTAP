clear all;

vector = [ -10.0000,    0.0000,    0.0000,    0.9886,    1.3974,    1.3974];

Kp_1=vector(1); % człon proporcjonalny
Ki_1=vector(2); % człon całkujący
Kd_1=vector(3);
Kp_2=vector(4); % człon proporcjonalny
Ki_2=vector(5); % człon całkujący
Kd_2=vector(6);

run("stale.m");

simulation_time = 2500; % Czas symulacji (sekundy)
T_p = 0.1; % Krok symulacji (sekundy)

tau_steps = tau / T_p;
tau_C_steps = tau_C / T_p;
k_min = max(tau_steps, tau_C_steps) + 1;
k_max = simulation_time / T_p + k_min;

% Trajektorie wejść procesu
T_H(1:k_max) = T_Hpp;
T_C(1:k_max) = T_Cpp;
T_D(1:k_max) = T_Dpp;
F_H(1:k_max) = F_Hpp;
F_Cin(1:k_max) = F_Cpp;
F_D(1:k_max) = F_Dpp;

T_zad(1:k_max) = 36.83;
T_zad(round(k_max/2):k_max) = 30.83 - 2;

h_zad(1:k_max) = 12.96;
h_zad(round(k_max/2):k_max) = 12.96 + 3;

% Stan i wyjścia procesu przed rozpoczęciem symulacji
x_1 = C * h_pp^3;
x_2 = C * h_pp^3 * T_pp;
x = [x_1, x_2]';
h(1:k_min) = h_pp;
T(1:k_min) = T_pp;
T_out(1:k_min) = T_pp;

% Wczytanie macierzy A, B, C, D
run("wspolczynniki_rownania_stanu.m");

% Symulacja (na podstawie równań stanu)
% Metoda RK2

e1(1:k_max) = 0;
e2(1:k_max) = 0;
e = 0;

% Kp_1 = 5;
% Ki_1 = 0.1;
% Kd_1 = 0.001;

r0_y1 = Kp_1+Ki_1*T_p+Kd_1/T_p;
r1_y1 = Kp_1+2*Kd_1/T_p;
r2_y1 = Kd_1/T_p;

% Kp_2 = 5;
% Ki_2 = 0.1;
% Kd_2 = 0.001;

r0_y2 = Kp_1+Ki_1*T_p+Kd_1/T_p;
r1_y2 = Kp_1+2*Kd_1/T_p;
r2_y2 = Kd_1/T_p;

for k = k_min:k_max

    e1(k) = T_zad(k) - T_out(k-1);
    F_Cin(k) = F_Cin(k-1) + r0_y1*e1(k) - r1_y1*e1(k-1) + r2_y1*e1(k-2);

    e2(k) = h_zad(k) - h(k-1);
    F_H(k) = F_H(k-1) + r0_y2*e2(k) - r1_y2*e2(k-1) + r2_y2*e2(k-2);

    e = e + abs(e1(k)) + abs(e2(k));
    
    [h, T_out, T] = obiekt(F_Cin, F_H, F_D, T_H, T_C, T_D, tau_C_steps, tau_steps, x, A, B, C, D, h, T, T_out, k, T_p);
end


figure(1)
hold on
plot(T_out(k_min:k_max))
plot(T_zad(k_min:k_max))
plot(h(k_min:k_max))
plot(h_zad(k_min:k_max))
hold off
legend(["T_{out}", "Tzad", "h", "hzad"])

figure(2)
hold on
plot(F_Cin(k_min:k_max))
plot(F_H(k_min:k_max))
hold off
legend(["Fcin", "Fh"])
