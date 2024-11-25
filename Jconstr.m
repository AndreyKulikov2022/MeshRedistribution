function [minus_detJ,ceq]=Jconstr(x,h)
[d_xi1,d_xi2] = diff_xi(x,h);
J=cat(2,d_xi1,d_xi2);
detJ=zeros(1,1,size(x,3),size(x,4));
for i=1:size(x,3)
    for j=1:size(x,4)
        detJ(:,:,i,j)=det(J(:,:,i,j));
    end
end
minus_detJ=-detJ;
ceq=0;
end