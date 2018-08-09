function M = feature_matrix_simulation( d1, d2, group )
% input
%   d1: dimention 1 of M
%   d2: dimention 2 of M
%   group

M = zeros( d1, d2);
for i = 1: d2
    M( :, i ) = column_simulation( d1, group );
end

end