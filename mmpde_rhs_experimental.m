function [result,bs,bders,as,aders,q_ali,B] = mmpde_rhs_experimental(x,metric,I,h,B)
%x=reshape(x,2,1,size(I,3),size(I,4));
tau=1;
[E_L,p,bs,bders,as,aders,q_ali,B]=variation_ali_experimental(x,metric,I, h,B);
result = E_L./p./tau;
% result(:,:,1,:)=0;
% result(:,:,end,:)=0;
% result(:,:,:,1)=0;
% result(:,:,:,end)=0;
% result=reshape(result,[],1);
end