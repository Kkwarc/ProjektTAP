clear all;
syms alpha C;
syms x_1 x_2 u_1 u_2 u_3
syms F_Dpp T_Hpp T_Cpp T_Dpp
syms x_1pp x_2pp;

%% Funkcje nieliniowe

f_1 = nthroot(x_1, 6);
f_2 = nthroot(x_1, 6) * x_2 / x_1;
f_3 = nthroot(x_1, 3);
f_4 = x_2 / x_1;

%% Pochodne cząstkowe

df_1dx_1 = diff(f_1, x_1);
df_2dx_1 = diff(f_2, x_1);
df_2dx_2 = diff(f_2, x_2);
df_3dx_1 = diff(f_3, x_1);
df_4dx_1 = diff(f_4, x_1);
df_4dx_2 = diff(f_4, x_2);

%% Funkcje zlinearyzowane

f_1L = subs(f_1, x_1, x_1pp) + subs(df_1dx_1, x_1, x_1pp) * (x_1 - x_1pp);
f_2L = subs(f_2, [x_1, x_2], [x_1pp, x_2pp]) + subs(df_2dx_1, [x_1, x_2],[ x_1pp, x_2pp]) * (x_1 - x_1pp) + subs(df_2dx_2, [x_1, x_2], [x_1pp, x_2pp]) * (x_2 - x_2pp);
f_3L = subs(f_3, x_1, x_1pp) + subs(df_3dx_1, x_1, x_1pp) * (x_1 - x_1pp);
f_4L = subs(f_4, [x_1, x_2], [x_1pp, x_2pp]) + subs(df_4dx_1, [x_1, x_2],[ x_1pp, x_2pp]) * (x_1 - x_1pp) + subs(df_4dx_2, [x_1, x_2], [x_1pp, x_2pp]) * (x_2 - x_2pp);


%% Równania stanu zlinearyzowane

dx_1 = u_1 + u_2 + F_Dpp - alpha/nthroot(C, 6) * f_1L;
dx_2 = T_Hpp * u_1 + T_Cpp * u_2 + F_Dpp * u_3 - alpha/nthroot(C, 6) * f_2L;
y_1 = 1 / nthroot(C, 3) * f_3L;
y_2 = f_4L;

%% Podstawienie punktu pracy i stałych
run("stale.m");
x_1pp = C * h_pp^3;
x_2pp = x_1pp * T_pp;

dx_1_subbed = subs(dx_1);
dx_2_subbed = subs(dx_2);
y_1_subbed = subs(y_1);
y_2_subbed = subs(y_2);