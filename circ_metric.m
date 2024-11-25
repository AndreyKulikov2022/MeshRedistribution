function M=circ_metric(x,v,R)
v_=[v(2);-v(1)];
ro=1+100*exp(-50*abs((x(1)-1/2)^2+(x(2)-1/2)^2-R^2));
M=ro^2 * v*(v')+v_*v_';
end