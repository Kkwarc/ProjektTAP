Tp = 10;% Krok symulacji (sekundy)
D = 100;

S = DMCstepmatrices(Tp, Tp*D);
[ny, nu, D] = size(S);

figure(10)
hold on
plot(reshape(S(1, 1, :), 1, D))
plot(reshape(S(1, 2, :), 1, D))
plot(reshape(S(2, 1, :), 1, D))
plot(reshape(S(2, 2, :), 1, D))
legend("s11", "s12", "s21", "s22")
title("Odpowiedzi skokowe")
hold off
print("odp_skokowe.eps","-depsc","-r400")
savefig("odp_skokowe.fig")
