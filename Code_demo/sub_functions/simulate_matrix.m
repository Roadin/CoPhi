function [ X, Y, G, E_org ] = simulate_matrix( m, n, p1, p2, group1, group2, dump_G, dump_E, noise)
%%
%E_org = X'GY
% E_org is m by n matrix
% X is m by p1 matrix
% Y is n by p2 matrix
% G is p1 by p2 matrix
%
% dump_G is percent of G we want to through away
% dump_E is percent of E we want to through away

% X = feature_matrix_simulation( m, p1-1, group1 );
% X = [ X sum( X, 2 )];
% Y = feature_matrix_simulation( n, p2-1, group2 );
% Y = [ Y sum( Y, 2 )];

 X = feature_matrix_simulation( m, p1, group1 );
Y =feature_matrix_simulation( n, p2, group2 );


G = normrnd( 0, 1, [ p1, p2 ]);
index = randsample( 1:1:p1*p2, floor( p1*p2*dump_G));
G( index ) = 0;

E_org = X*G*Y' + normrnd( 0, noise, [ m, n]);
index = randsample( 1:1:m*n, floor( m*n*dump_E));
E_org( index ) = nan;

end
