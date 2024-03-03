% constans
C = 0.6;
alpha = 15;

% working point
tauC = 90;
tau = 40;

delay = max(tauC, tau);
sim_end = 1000+delay;


Tc(1:1:sim_end) = 19;

Th(1:1:sim_end) = 73;

Td(1:1:sim_end) = 24;

Fc(1:1:sim_end) = 28;

Fh(1:1:sim_end) = 17;

Fd(1:1:sim_end) = 9;

h(1:1:sim_end) = 0;

T(1:1:sim_end) = 0;

V(1:1:sim_end) = 0;
F(1:1:sim_end) = 0;