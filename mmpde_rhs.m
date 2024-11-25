function [result,bs,bders,as,aders,q_ali] = mmpde_rhs(x,metric,I,h)
%x=reshape(x,2,1,size(I,3),size(I,4));
tau=1;
[E_L,p,bs,bders,as,aders,q_ali]=variation_ali(x,metric,I, h);
result = E_L./p./tau;
 result(:,:,1,:)=0;
 result(:,:,end,:)=0;
 result(:,:,:,1)=0;
 result(:,:,:,end)=0;
% result=reshape(result,[],1);
end