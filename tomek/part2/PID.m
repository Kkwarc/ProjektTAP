function e = PID(vector)
Kp_1=vector(1); % człon proporcjonalny
Ki_1=vector(2); % człon całkujący
Kd_1=vector(3);
Kp_2=vector(4); % człon proporcjonalny
Ki_2=vector(5); % człon całkujący
Kd_2=vector(6);

% Stałe
C = 0.6;
alpha = 15;
tau_C = 90;
tau = 40;

% Punkt pracy - wejścia procesu
T_Cpp = 19;
T_Hpp = 73;
T_Dpp = 24;
F_Cpp = 28;
F_Hpp = 17;
F_Dpp = 9;

% Punkt pracy - wyjścia procesu
h_pp = 12.96;
T_pp = 36.83;

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
T_zad(round(k_max/3):k_max) = 36.83 - 2;
T_zad(round(2*k_max/3):k_max) = 36.83 + 2;

h_zad(1:k_max) = 12.96;
h_zad(round(k_max/3):k_max) = 12.96 + 1;
h_zad(round(2*k_max/3):k_max) = 12.96 - 1;

% Stan i wyjścia procesu przed rozpoczęciem symulacji
F(1:k_min) = alpha * sqrt(h_pp);
F_C(1:k_min) = F_Cpp;
h(1:k_min) = h_pp;
V(1:k_min) = C * h_pp^3;
T(1:k_min) = T_pp;
VT(1:k_min) = V(1) * T(1);
T_out(1:k_min) = T_pp;

e1(1:k_max) = 0;
e2(1:k_max) = 0;
e = 0;

r0_y1 = Kp_1+Ki_1*T_p+Kd_1/T_p;
r1_y1 = Kp_1+2*Kd_1/T_p;
r2_y1 = Kd_1/T_p;

r0_y2 = Kp_2+Ki_2*T_p+Kd_2/T_p;
r1_y2 = Kp_2+2*Kd_2/T_p;
r2_y2 = Kd_2/T_p;

for k = k_min:k_max

    e1(k) = T_zad(k) - T_out(k-1);
    F_Cin(k) = F_Cin(k-1) + r0_y1*e1(k) - r1_y1*e1(k-1) + r2_y1*e1(k-2);
    if F_Cin(k) < 0
        F_Cin(k) = 0;
    end

    e2(k) = h_zad(k) - h(k-1);
    F_H(k) = F_H(k-1) + r0_y2*e2(k) - r1_y2*e2(k-1) + r2_y2*e2(k-2);
    if F_H(k) < 0
        F_H(k) = 0;
    end

    e = e + abs(e1(k)) + abs(e2(k));
    
    [F_C, V, VT, T, F, h, T_out] = obiekt(F_Cin, F_H, F_D, F_C, T_H, T_C, T_D, T_out, h, C, alpha, tau_C_steps, tau_steps, V, VT, T, F, T_p, k);
end

disp(e)
disp(vector)

end