function F=functional_shifts(x, metric, Q, K, delta)
F = functional(x, metric, Q, K);
F = ones(size(x))*F;
for i=2:size(x,3)-1
    for j = 2:size(x,4)-1
    x_delta = x;
    y_delta = x;
    x_delta(1,1,i,j) = x(1,1,i,j)+delta;
    y_delta(2,1,i,j) = x(2,1,i,j)+delta;
    F(1,1,i,j) = functional(x_delta,metric, Q, K);
    F(2,1,i,j) = functional(y_delta,metric, Q, K);
    end
end

end