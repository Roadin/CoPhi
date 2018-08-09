function [v] = soft_thresh(alpha, beta)
% apply soft-threshold with given alpha and beta

% alpha and beta can be vector or matrix

v = zeros(size(alpha));
alpha_abs = abs(alpha);
nonZero = alpha_abs > beta;
v(nonZero) = sign(alpha(nonZero)) .* (alpha_abs(nonZero) - beta(nonZero));

end

