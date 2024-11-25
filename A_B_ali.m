function [As, Bs] = A_B_ali(x,metric,I, h)
%% M, I - 4-D arrays, first two are for the matrix itself, third and second is for 2D domain 
% x - coodrdinates
% M - metric tensor, I - identity matrices.
% h - uniform grid step
[d_xi1,d_xi2] = diff_xi(x,h);
J=cat(2,d_xi1,d_xi2);
Ji=pageinv(J);
M=zeros(2,2,size(x,3),size(x,4));
ro=zeros(1,1,size(x,3),size(x,4));
for i=1:size(x,3)
    for j=1:size(x,4)
        M(:,:,i,j)=metric(x(:,:,i,j));
        ro(:,:,i,j)=sqrt(det(M(:,:,i,j)));
    end
end
Mi=pageinv(M);
[rod_xi1,rod_xi2] = diff_xi(ro,h);
rods={rod_xi1,rod_xi2};
[Mid_xi1,Mid_xi2] = diff_xi(Mi,h);
Mids={Mid_xi1,Mid_xi2};
As=cell(2,2);
Bs={zeros(size(I)),zeros(size(I))};
for i=1:2
    for j=1:2
        aMa=pagemtimes(Ji(i,:,:,:), ...
            pagemtimes(Mi,'none',Ji(j,:,:,:),'transpose'));
       
        As{i,j}=pagemtimes(I,pagemtimes(ro,aMa));

        aMida=pagemtimes(Ji(i,:,:,:), ...
            pagemtimes(Mids{i},'none',Ji(j,:,:,:),'transpose'));
        Bs{j}=Bs{j}-pagemtimes(I, ...
            pagemtimes(aMa,rods{i})+pagemtimes(ro,aMida));
    end
end
end