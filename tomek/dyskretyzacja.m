clear all;

run("wspolczynniki_rownania_stanu.m");

ss_model = ss(A, B, C, D);

metoda = "zoh"; % Metoda dyskretyzacji
T_p = 15.0; % Czas próbkowania

ss_model_dyskretny = c2d(ss_model, T_p, metoda);

% Output poniższych poleceń można przekopiować wygodnie do latexa i są od
% razu ładnie sformatowane macierze
disp("\[\begin{aligned}")
disp("A &=")
disp(replace(replace(latex(vpa(ss_model_dyskretny.A, 4)), "(", "["), ")", "]"))
disp("\quad")
disp("B =")
disp(replace(replace(latex(vpa(ss_model_dyskretny.B, 4)), "(", "["), ")", "]"))
disp("\\")
disp("C &=")
disp(replace(replace(latex(vpa(ss_model_dyskretny.C, 4)), "(", "["), ")", "]"))
disp("\quad")
disp("D =")
disp(replace(replace(latex(vpa(ss_model_dyskretny.D, 4)), "(", "["), ")", "]"))
disp("\end{aligned}\]")

A = ss_model_dyskretny.A;
B = ss_model_dyskretny.B;
C = ss_model_dyskretny.C;
D = ss_model_dyskretny.D;

save wspolczynniki_dyskretnego_rownania_stanu A B C D T_p