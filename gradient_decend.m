function x_new = gradient_decend(x, metric, Q, K)

delta = 0.1;
tau = 0.003;

F_plus = functional_shifts(x, metric, Q, K, delta); 
F_minus = functional_shifts(x, metric, Q, K, -delta);

%use hessian pre-conditioning?
dF = (F_plus - F_minus)/(2*delta);

x_new = x - tau*dF;

end

