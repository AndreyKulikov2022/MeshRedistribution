function M=sin_metric(x,v,N)
v_=[v(2);-v(1)];
ro=1+100*exp(-50*abs(x(2)/N-1/2-1/4*sin(2*pi*x(1)/N)));
M=ro^2 * v*(v')+v_*v_';
end