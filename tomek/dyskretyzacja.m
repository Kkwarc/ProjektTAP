clear all;

%% Równania stanu

run("wspolczynniki_rownania_stanu.m");

ss_model = ss(A, B, C, D);

metoda = "zoh"; % Metoda dyskretyzacji
T_p = 10.0; % Czas próbkowania

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

%% Transmitancja

tf_model = tf(ss_model_dyskretny);
syms z;

% Heavy wizardry begins here
disp("\[\begin{aligned}")
for i = 1:2
    for j = 1:3
        disp("G_{" + i + "," + j + "}(z) &=")
        num = tf_model.Numerator(i, j);
        num = num{1, 1};
        den = tf_model.Denominator(i, j);
        den = den{1, 1};
        TF=tf(num,den,1,'variable','z');
        fracparts=regexp(evalc('TF'),'([^\n]*)\n[ ]*-[-]+[ ]*\n([^\n]*)','tokens');
        TFlatex=['\frac{' fracparts{1}{1} '}{' fracparts{1}{2} '}'];
        disp(TFlatex);
        disp("\\ \\")
    end
end
disp("K_{1} &=")
num = tf_model.Numerator(1, 4);
num = num{1, 1};
den = tf_model.Denominator(1, 4);
den = den{1, 1};
disp(latex(vpa(num(1)/den(1),4)))
disp("\\ \\")
disp("K_{2} &=")
num = tf_model.Numerator(2, 4);
num = num{1, 1};
den = tf_model.Denominator(2, 4);
den = den{1, 1};
disp(latex(vpa(num(1)/den(1),4)))
disp("\end{aligned}\]")
% Heavy wizardry end.

save transmitancja_dyskretna tf_model
