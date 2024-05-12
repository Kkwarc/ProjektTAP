% function [E, h1, h2, h2_zad, F1, Fd]=DMC(wektor, zaklocenia, rysowanie)

clear all
close all

Tp = 10;% Krok symulacji (sekundy)
D = 100;

start = D+1;
simulation_time = 10000; % Czas symulacji (sekundy)
poczatek = start; %chwila k w której zmienia sie wartość zadana
N = 35;
Nu = 20;
lambda = [200, 20];
phi = [1, 10];

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

T_zad(1:simulation_time) = T_pp;
T_zad(start:start+round((simulation_time-start)/3)) = T_pp;
T_zad(round(start+(simulation_time-start)/3):start+round(2*(simulation_time-start)/3)) = T_pp - 5;
T_zad(start+round(2*(simulation_time-start)/3):simulation_time) = T_pp + 2;

h_zad(1:simulation_time) = h_pp;
h_zad(start:start+round((simulation_time-start)/3)) = h_pp;
h_zad(round(start+(simulation_time-start)/3):start+round(2*(simulation_time-start)/3)) = h_pp + 4;
h_zad(start+round(2*(simulation_time-start)/3):simulation_time) = h_pp - 1;

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

% figure(10)
% hold on
% plot(reshape(S(1, 1, :), 1, D))
% plot(reshape(S(1, 2, :), 1, D))
% plot(reshape(S(2, 1, :), 1, D))
% plot(reshape(S(2, 2, :), 1, D))
% % plot(reshape(S(3, 1, :), 1, D))
% % plot(reshape(S(3, 2, :), 1, D))
% % plot(reshape(S(4, 1, :), 1, D))
% % plot(reshape(S(4, 2, :), 1, D))
% legend("s11", "s11", "s21", "s22")
% hold off

lambda_mat = [lambda(1), 0; 0, lambda(2)];
LAMBDA = kron(eye(Nu), lambda_mat);

phi_mat = [phi(1), 0; 0, phi(2)];
PHI = kron(eye(N), phi_mat);

[M, MP] = DMCmatrices(S, N, Nu);
% Sz = DMCstepmatricesZ(Tp, timefinal);
% MzP = DMCmatrixMzP(Sz, N);

K = (M'*PHI*M+LAMBDA)^(-1)*M'*PHI;

DU_p = zeros((D-1)*nu, 1);
for k=start:simulation_time
%     disp(k)
    %symulacja obiektu
    [F_C, V, VT, T, F, h, T_out] = obiekt(F_Cin, F_H, F_D, F_C, T_H, T_C, T_D, T_out, h, C, alpha, tau_C_steps, tau_steps, V, VT, T, F, Tp, k);

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

    F_Cin(k) = F_Cin(k-1) + DU(1);

    if F_Cin(k) > F_Cpp*1.4
        F_Cin(k) = F_Cpp*1.4;
    elseif F_Cin(k) < F_Cpp*0.6
        F_Cin(k) = F_Cpp*0.6;
    end

    F_H(k) = F_H(k-1) + DU(2);
    
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
stairs(h(1:end))
plot(h_zad(1:end))
stairs(T_out(1:end))
plot(T_zad(1:end))
title("Wyjście")
legend("Wyjście h", "h zadana", "Wyjscie Tout", "t zadana")
xlabel("chwila k")
ylabel("wartość wyjścia")
hold off

subplot(2,1,2)
hold on
stairs(F_Cin(1:end))
stairs(F_H(1:end))
legend("sterowanie Fcin", "sterowanie Fh")
xlabel("chwila k")
ylabel("wartość sterowania")
title("Sterowanie")
print("DMC_single.eps","-depsc","-r400")
