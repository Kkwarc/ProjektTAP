clear all;
close all;

run('stale.m');

h0 = 12.96;
T0 = 36.83;

for t=delay:sim_end
    F(t) = alpha*sqrt(h0) + alpha*(1 / (2 * sqrt(h0))) * (h(t-1) - h0);
    V(t) = V(t-1) + Fh(t-1) + Fc(t-1) + Fd(t-1) - F(t);

    T(t) = T(t-1) + (Fh(t-1)*Th(t-1)+Fc(t-1)*Tc(t-1)+Fd(t-1)*Td(t-1)-F(t-1)*T(t-1))/V(t);

    h(t) = h0 + (V(t)-C*h0^3)/(3*C*h0^2);
    Tout(t) = T(t-tau);
end

figure(1)
hold on
plot(h)
plot(Tout)
legend("h", "Tout")