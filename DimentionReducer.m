% For better representation of data, we need to convert 2D arrays to 1D
% here...
function [one_d_arr]=DimentionReducer(two_d_arr)
sz = size(two_d_arr);
one_d_arr = zeros(1,sz(1) * sz(2));
index = 1;
for i=1:sz(1)
    for j=1:sz(2)
        one_d_arr(index)=two_d_arr(i,j);
        index = index + 1;
    end
end
end
