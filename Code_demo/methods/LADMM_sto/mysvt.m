function E = mysvt(Y, x)
%% singular value thresholding operator for matrix Y by thretholding parameter x
% sv = 5;
% [m, n]=size(Y);
% opt.tol = tol2;%precision for computing the partial SVD
% opt.p0 = ones(n,1);
% [U, S, V] = lansvd(Y, m, n, sv, 'L', opt);
%     %[U, S, V] = lansvd(M, n, n, sv, 'L');
%     %[U, S, V] = svd(M,'econ');
%       
%     S = diag(S);
%     svp = length(find(S>x));
%     if svp < sv
%         sv = min(svp + 1, n);
%     else
%         sv = min(svp + round(0.05*n), n);
%     end
%     
%     if svp>=1
%         S = S(1:svp) - x;
%     else
%         svp = 1;
%         S = 0;
%     end
% 
%     A.U = U(:, 1:svp);
%     A.s = S;
%     A.V = V(:, 1:svp);
%     
%     E = A.U*diag(A.s)*A.V';
    

[S, V, D] = svd(Y);
v = diag(V);
[V_row, V_col] = size(V);
x = x * ones(size(v));
v_new = zeros(size(v));
nonZero = v > x;
v_new(nonZero) = v(nonZero) - x(nonZero);
if V_row<V_col
E = S * [diag(v_new), zeros(V_row, V_col-V_row)] * D';
else E = S * [diag(v_new); zeros(V_row-V_col, V_col)] * D';
end