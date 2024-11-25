function x=non_conservative(x,metric,I,h,dt,n_steps)
for i=1:n_steps
    rhs=mmpde_rhs(x,metric,I,h);
x(:,:,2:end-1,2:end-1)=x(:,:,2:end-1,2:end-1)+dt*rhs(:,:,2:end-1,2:end-1);
end