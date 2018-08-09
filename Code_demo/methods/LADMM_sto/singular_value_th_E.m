function E_new = singular_value_th_E(F, E, G_new, C_new, Omega, X, Y, lambda_E, multiplier_1, multiplier_2, tau, beta)
%% solve E via SVT
gradient_E_1 = Omega .* (E - F + multiplier_1 / beta);
gradient_E_2 = E - C_new - X' * G_new * Y + multiplier_2 / beta;
Z = E - 0.5 / tau * (gradient_E_1 + gradient_E_2);

% E_new = mysvt(Z, lambda_E / (2 * tau *  beta));

[E_new, ~] = Pro2TraceNorm(Z, 2*tau*beta);
end