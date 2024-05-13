% function [E, h1, h2, h2_zad, F1, Fd]=DMC(wektor, zaklocenia, rysowanie)

clear all
close all

Tp = 10;% Krok symulacji (sekundy)
D = 69;

start = D+1;
simulation_time = 10000; % Czas symulacji (sekundy)
poczatek = start; %chwila k w której zmienia sie wartość zadana
N = 25;
Nu = 5;
lambda = [1, 1];
phi = [1, 1];

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

tau_steps = tau / Tp;
tau_C_steps = tau_C / Tp;

% Trajektorie wejść procesu
T_H(1:simulation_time) = T_Hpp;
T_C(1:simulation_time) = T_Cpp;
T_D(1:simulation_time) = T_Dpp;
F_H(1:simulation_time) = F_Hpp;
F_Cin(1:simulation_time) = F_Cpp;
F_D(1:simulation_time) = F_Dpp;

minimal_level = 0.9;
maximal_level = 1.1;

F_D(1:simulation_time) = F_Dpp;
F_D(start:start+round((5/9)*(simulation_time-start))) = F_Dpp;
F_D(round((5/9)*(simulation_time-start)):round((8/9)*(simulation_time-start))) = F_Dpp * minimal_level;
F_D(round((8/9)*(simulation_time-start)):simulation_time) = F_Dpp * maximal_level;

T_zad(1:simulation_time) = T_pp;
T_zad(start:start+round((3/9)*(simulation_time-start))) = T_pp;
T_zad(round((3/9)*(simulation_time-start)):round((6/9)*(simulation_time-start))) = T_pp * minimal_level;
T_zad(round((6/9)*(simulation_time-start)):simulation_time) = T_pp * maximal_level;

h_zad(1:simulation_time) = h_pp;
h_zad(start:round((4/9)*(simulation_time-start))) = h_pp;
h_zad(round((4/9)*(simulation_time-start)):round((7/9)*(simulation_time-start))) = h_pp * minimal_level;
h_zad(round((7/9)*(simulation_time-start)):simulation_time) = h_pp * maximal_level;

% Stan i wyjścia procesu przed rozpoczęciem symulacji
F(1:simulation_time) = alpha * sqrt(h_pp);
F_C(1:simulation_time) = F_Cpp;
h(1:simulation_time) = h_pp;
V(1:simulation_time) = C * h_pp^3;
T(1:simulation_time) = T_pp;
VT(1:simulation_time) = V(1) * T(1);
T_out(1:simulation_time) = T_pp;

e = 0;

S = DMCstepmatrices(Tp, Tp*D);
[ny, nu, D] = size(S);

lambda_mat = [lambda(1), 0; 0, lambda(2)];
LAMBDA = kron(eye(Nu), lambda_mat);

phi_mat = [phi(1), 0; 0, phi(2)];
PHI = kron(eye(N), phi_mat);

[M, MP] = DMCmatrices(S, N, Nu);

K = (M'*PHI*M+LAMBDA)^(-1)*M'*PHI;

DU_p = zeros((D-1)*nu, 1);
for k=start:simulation_time
%     disp(k)
    %symulacja obiektu
    [F_C, V, VT, T, F, h, T_out] = obiekt(F_Cin, F_H, F_D, F_C, T_H, T_C, T_D, T_out, h, C, alpha, tau_C_steps, tau_steps, V, VT, T, F, Tp, k);

    %Obliczenie DU_p
    for d=1:(D-1)
        DU_p(2*d-1) = F_H(k-d) - F_H(k-d-1);
        DU_p(2*d) = F_Cin(k-d) - F_Cin(k-d-1);
    end

    %Pomiar wyjścia
    Y = ones(N*ny, 1);
    Y(1:2:end) = Y(1:2:end) * h(k);
    Y(2:2:end) = Y(2:2:end) * T_out(k);

    Y_zad = ones(N*ny, 1);
    Y_zad(1:2:end) = Y_zad(1:2:end) * h_zad(k);
    Y_zad(2:2:end) = Y_zad(2:2:end) * T_zad(k);

    %Obliczenie sterowania
    DU = K*(Y_zad-Y-MP*DU_p);

    F_Cin(k) = F_Cin(k-1) + DU(2);

    if F_Cin(k) > F_Cpp*1.4
        F_Cin(k) = F_Cpp*1.4;
    elseif F_Cin(k) < F_Cpp*0.6
        F_Cin(k) = F_Cpp*0.6;
    end

    F_H(k) = F_H(k-1) + DU(1);
    
    if F_H(k) > F_Hpp*1.4
        F_H(k) = F_Hpp*1.4;
    elseif F_H(k) < F_Hpp*0.6
        F_H(k) = F_Hpp*0.6;
    end
    e = e + (T_zad(k)-T_out(k))^2 + (h_zad(k)-h(k))^2; 
end
disp(e)

figure(2)
subplot(2,1,1)
hold on
stairs(h(start:end))
plot(h_zad(start:end))
stairs(T_out(start:end))
plot(T_zad(start:end))
plot(F_D(start:end))
title("Wyjście")
legend("Wyjście h", "h zadana", "Wyjscie Tout", "t zadana", "zak")
xlabel("chwila k")
ylabel("wartość wyjścia")
hold off

subplot(2,1,2)
hold on
stairs(F_Cin(start:end))
stairs(F_H(start:end))
legend("sterowanie Fcin", "sterowanie Fh")
xlabel("chwila k")
ylabel("wartość sterowania")
title("Sterowanie")
print("DMC_analityczny.eps","-depsc","-r400")
