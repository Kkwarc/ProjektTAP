% function [E, h1, h2, h2_zad, F1, Fd]=DMC(wektor, zaklocenia, rysowanie)

clear all
close all

Tp = 0.1;% Krok symulacji (sekundy)
D = 5000;

start = D+1;
simulation_time = 10000; % Czas symulacji (sekundy)
poczatek = start; %chwila k w której zmienia sie wartość zadana
N = 120;
Nu = 20;
lambda = 10;

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
k_min = max(tau_steps, tau_C_steps) + 1;
k_max = simulation_time / Tp + k_min;

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
F(1:k_max) = alpha * sqrt(h_pp);
F_C(1:k_max) = F_Cpp;
h(1:k_max) = h_pp;
V(1:k_max) = C * h_pp^3;
T(1:k_max) = T_pp;
VT(1:k_max) = V(1) * T(1);
T_out(1:k_max) = T_pp;

e = 0;

S = DMCstepmatrices(Tp, Tp*D);
[ny, nu, D] = size(S);

figure(1)
hold on
plot(reshape(S(1, 1, :), 1, D))
plot(reshape(S(1, 2, :), 1, D))
plot(reshape(S(2, 1, :), 1, D))
plot(reshape(S(2, 2, :), 1, D))
% plot(reshape(S(3, 1, :), 1, D))
% plot(reshape(S(3, 2, :), 1, D))
% plot(reshape(S(4, 1, :), 1, D))
% plot(reshape(S(4, 2, :), 1, D))
legend("s11", "s11", "s21", "s22")
hold off

L = eye(Nu*nu, Nu*nu)*lambda;

[M, MP] = DMCmatrices(S, N, Nu);
% Sz = DMCstepmatricesZ(Tp, timefinal);
% MzP = DMCmatrixMzP(Sz, N);

K = (M'*M+L)^(-1)*M';

DU_p = zeros((D-1)*nu, 1);
for k=start:simulation_time
    disp(k)
    %symulacja obiektu
    [F_C, V, VT, T, F, h, T_out] = obiekt(F_Cin, F_H, F_D, F_C, T_H, T_C, T_D, T_out, h, C, alpha, tau_C_steps, tau_steps, V, VT, T, F, Tp, k);
    
    if T_out(k) > 100
        T_out(k) = 100;
    elseif T_out(k) < 5
        T_out(k) = 5;
    end

    %Obliczenie DU_p
    for d=1:(D-1)
        DU_p(2*d-1) = F_Cin(k-d) - F_Cin(k-d-1);
        DU_p(2*d) = F_H(k-d) - F_H(k-d-1);
    end

    %Pomiar wyjścia
    Y = ones(N*ny, 1);
    Y(1:2:end) = Y(1:2:end) * T_out(k);
    Y(2:2:end) = Y(2:2:end) * h(k);

    Y_zad = ones(N*ny, 1);
    Y_zad(1:2:end) = Y_zad(1:2:end) * T_zad(k);
    Y_zad(2:2:end) = Y_zad(2:2:end) * h_zad(k);

    %Obliczenie sterowania
    DU = K*(Y_zad-Y-MP*DU_p);

    F_Cin(k) = DU(1);

    if F_Cin > 100
        F_Cin = 1;
    elseif F_Cin < 0
        F_Cin = 0;
    end

    F_H(k) = DU(3);

    if F_H > 100
        F_H = 1;
    elseif F_H < 0
        F_H = 0;
    end
    e = e + (T_zad(k)-T_out(k))^2 + (h_zad(k)-h(k))^2; 
end

subplot(2,1,1)
hold on
stairs(h(start:end))
plot(h_zad(start:end))
stairs(T_out(start:end))
plot(T_zad(start:end))
title("Wyjście")
legend("Wyjście h", "h zadana", "Wyjscie Tout", "t zadana")
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
print("DMC_single.eps","-depsc","-r400")
