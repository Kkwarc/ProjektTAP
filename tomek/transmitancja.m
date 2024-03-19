clear all;

% Wczytanie macierzy A, B, C, D
run("wspolczynniki_rownania_stanu.m");

syms s;

G = C * inv(s * eye(2) - A)*B + D;

for i = 1:size(G, 1)
    for j = 1:size(G, 2)
        [N, D] = numden(G(i, j));
        N = collect(N, s);
        D = collect(D, s);
        Dc = coeffs(D, s);
        NN = N/Dc(1);
        ND = D/Dc(1);
        disp(latex(vpa(NN/ND)))
    end

end