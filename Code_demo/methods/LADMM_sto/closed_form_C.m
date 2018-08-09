function C_new = closed_form_C(beta, E, G, X, Y, multiplier_2)
C_new = (1/(beta+1))*(beta*( E-X'* G * Y)+ multiplier_2);
end