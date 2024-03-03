clear all;
close all;

run('stale.m');

for t=delay:sim_end
    F(t) = alpha * sqrt(h(t-1));
    V(t) = V(t-1) + Fh(t-1) + Fc(t-1) + Fd(t-1) - F(t);

    T(t) = T(t-1) + (Fh(t-1)*Th(t-1)+Fc(t-1)*Tc(t-1)+Fd(t-1)*Td(t-1)-F(t-1)*T(t-1))/V(t);

    h(t) = nthroot(V(t)/C, 3);
    Tout(t) = T(t-tau);
end

figure(1)
hold on
plot(h)
plot(Tout)
legend("h", "Tout")