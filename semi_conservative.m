function x=semi_conservative(x,metric,I,h,dt,n_steps)
[d_xi1,d_xi2] = diff_xi(x,h);
J=cat(2,d_xi1,d_xi2);
detJ=zeros(1,1,size(x,3),size(x,4));
for i=1:size(x,3)
    for j=1:size(x,4)
        detJ(:,:,i,j)=det(J(:,:,i,j));
    end
end
for i=1:n_steps
    [rhs, rhsJ]=semi_cons_rhs(x,metric,detJ,h);
    x=x+dt*rhs;
    detJ=detJ+dt*rhsJ;
end