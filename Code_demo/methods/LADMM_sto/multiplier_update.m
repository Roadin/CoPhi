function [M_1_new, M_2_new] = multiplier_update(M_1, M_2, E, F, G, X, Y, C, Omega, beta)
M_1_new = M_1 + beta * (Omega .*(E - F));
M_2_new = M_2 + beta * (E - X' * G * Y - C);
end