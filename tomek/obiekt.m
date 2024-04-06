function [F_C, V, VT, T, F, h, T_out] = obiekt(F_Cin, F_H, F_D, F_C, T_H, T_C, T_D, T_out, h, C, alpha, tau_C_steps, tau_steps, V, VT, T, F, T_p, k)
    F_C(k) = F_Cin(k - tau_C_steps);
    dVdt_1 = F_H(k-1) + F_C(k-1) + F_D(k-1) - F(k-1);
    dVTdt_1 = F_H(k-1) * T_H(k-1) + F_C(k-1) * T_C(k-1) + F_D(k-1) * T_D(k-1) - F(k-1) * T(k-1);
    V_temp = V(k-1) + 0.5 * T_p * dVdt_1;
    VT_temp = VT(k-1) + 0.5 * T_p * dVTdt_1;
    T_temp = VT_temp / V_temp;
    h_temp = nthroot(V_temp / C, 3);
    F_temp = alpha * sqrt(h_temp);
    dVdt_2 = F_H(k-1) + F_C(k-1) + F_D(k-1) - F_temp;
    dVTdt_2 = F_H(k-1) * T_H(k-1) + F_C(k-1) * T_C(k-1) + F_D(k-1) * T_D(k-1) - F_temp * T_temp;
    V(k) = V(k-1) + T_p * dVdt_2;
    VT(k) = VT(k-1) + T_p * dVTdt_2;
    T(k) = VT(k) / V(k);
    F(k) = alpha * sqrt(h(k-1));
    h(k) = nthroot(V(k-1) / C, 3); 
    T_out(k) = T(k - tau_steps);
end