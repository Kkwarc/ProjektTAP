function [h, T_out, T] = obiekt(F_Cin, F_H, F_D, T_H, T_C, T_D, tau_C_steps, tau_steps, x, A, B, C, D, h, T, T_out, k, T_p)
F_C = F_Cin(k-tau_C_steps);
u_1 = F_H(k);
u_2 = F_C;
u_3 = F_D(k);
u_4 = F_H(k) * T_H(k);
u_5 = F_C * T_C(k);
u_6 = F_D(k) * T_D(k);
u = [u_1, u_2, u_3, u_4, u_5, u_6, 1]';

k_1 = A*x + B*u;
x_temp = x + 0.5 * T_p * k_1;
k_2 = A*x_temp + B*u;
x = x + T_p * k_2;
y = C*x + D*u;

h(k) = y(1);
T(k) = y(2);
T_out(k) = T(k-tau_steps);
end