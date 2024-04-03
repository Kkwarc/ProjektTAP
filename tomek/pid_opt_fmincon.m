%Ograniczenia dolne parametrów
lb = [-10, 0.00000001, 0.0000001, -10, 0.0000001, 0.0000001];
% %Ograniczenia górne parametrów
ub = [10, 10, 10, 10, 10, 10];
%nastawy startowe
nastawy_startowe_PID = [1, 1, 1, 1, 1, 1];

disp(length(nastawy_startowe_PID))
optymalne_nastawy_PID = fmincon(@PID,nastawy_startowe_PID, [], [], [], [], lb, ub);
disp(optymalne_nastawy_PID)