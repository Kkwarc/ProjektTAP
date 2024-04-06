clear all;
close all;
clc;

%Ograniczenia dolne parametrów
lb = [-2, -2, 0, -2, -2, 0];
% %Ograniczenia górne parametrów
ub = [2, 2, 2, 2, 2, 2];
%nastawy startowe
nastawy_startowe_PID = [-1, 0, 0, 1, 0, 0];

disp(length(nastawy_startowe_PID))
optymalne_nastawy_PID = fmincon(@PID,nastawy_startowe_PID, [], [], [], [], lb, ub);
disp(optymalne_nastawy_PID)