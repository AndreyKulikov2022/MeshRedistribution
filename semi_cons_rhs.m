function [result, dJ]=semi_cons_rhs(x,metric,detJ,h)
result=zeros(size(x));
dJ=zeros(size(detJ));
[d_xi1,d_xi2] = diff_xi(x,h);
J=cat(2,d_xi1,d_xi2);
% detJ=zeros(1,1,size(x,3),size(x,4));
% for i=1:size(x,3)
%     for j=1:size(x,4)
%         detJ(:,:,i,j)=det(J(:,:,i,j));
%     end
% end
Ja1=cat(1,d_xi2(2,:,:,:),-d_xi2(1,:,:,:));
Ja2=cat(1,-d_xi1(2,:,:,:),d_xi1(1,:,:,:));
M=zeros(2,2,size(x,3),size(x,4));
ro=zeros(1,1,size(x,3),size(x,4));
for i=1:size(x,3)
    for j=1:size(x,4)
        M(:,:,i,j)=metric(x(:,:,i,j));
        ro(:,:,i,j)=sqrt(det(M(:,:,i,j)));
    end
end
Mi=pageinv(M);
roMidJ=ro.*Mi./detJ;
Ja1roma1=pagemtimes(pagemtimes(Ja1,'transpose',roMidJ,'none'),Ja1);
Ja1roma2=pagemtimes(pagemtimes(Ja1,'transpose',roMidJ,'none'),Ja2);
Ja2roma2=pagemtimes(pagemtimes(Ja2,'transpose',roMidJ,'none'),Ja2);

[d_xi1a1a1,~]=diff_xi(Ja1roma1,h);
[d_xi1a2a1,d_xi2a2a1]=diff_xi(Ja1roma2,h);
[~,d_xi2a2a2]=diff_xi(Ja2roma2,h);

Bs{1}=-(d_xi1a1a1+d_xi2a2a1)./detJ;
Bs{2}=-(d_xi1a2a1+d_xi2a2a2)./detJ;
Js=cell(2,1);
d_xi1_fb = zeros(size(x));
d_xi1_forward=(x(:,:,2:end-1,3:end)-x(:,:,2:end-1,2:end-1))/(h);
d_xi1_backward=(x(:,:,2:end-1,2:end-1)-x(:,:,2:end-1,1:end-2))/(h);
ids1=Bs{1}(1,1,:,:)>0;
d_xi1_fb(:,:,2:end-1,2:end-1)=d_xi1_forward.*ids1(:,:,2:end-1,2:end-1)+d_xi1_backward.*(~ids1(:,:,2:end-1,2:end-1));
Js{1}=d_xi1_fb;
d_xi2_fb = zeros(size(x));
d_xi2_forward=(x(:,:,3:end,2:end-1)-x(:,:,2:end-1,2:end-1))/(h);
d_xi2_backward=(x(:,:,2:end-1,2:end-1)-x(:,:,1:end-2,2:end-1))/(h);
ids2=Bs{2}(1,1,:,:)>0;
d_xi2_fb(:,:,2:end-1,2:end-1)=d_xi2_forward.*ids2(:,:,2:end-1,2:end-1)+d_xi2_backward.*(~ids2(:,:,2:end-1,2:end-1));
Js{2}=d_xi2_fb;

% Balancing function
frob2=@(F) sum(sum(F.*F,1),2);
p=1;%sqrt(frob2(d_xi1a1a1)+frob2(d_xi1a2a1)+frob2(d_xi2a2a1)+frob2(d_xi2a2a2)).*detJ;

result=-(Js{1}.*(d_xi1a1a1+d_xi2a2a1)+...
    Js{2}.*(d_xi1a2a1+d_xi2a2a2))./detJ./p;

result(:,:,1,:)=0;
result(:,:,end,:)=0;
result(:,:,:,1)=0;
result(:,:,:,end)=0;
% 
[cd1,cd2]=diff_central(x,h);
Jxit=-cd2(2,:,:,:).*(result(1,:,2:end-1,1:end-1)+result(1,:,2:end-1,2:end))/2+...
cd2(1,:,:,:).*(result(2,:,2:end-1,1:end-1)+result(2,:,2:end-1,2:end))/2;
Jetat=cd1(2,:,:,:).*(result(1,:,1:end-1,2:end-1)+result(1,:,2:end,2:end-1))/2-...
cd1(1,:,:,:).*(result(2,:,1:end-1,2:end-1)+result(2,:,2:end,2:end-1))/2;
dJ(:,:,2:end-1,2:end-1)=-(Jxit(:,:,:,2:end)-Jxit(:,:,:,1:end-1))/h-(Jetat(:,:,2:end,:)-Jetat(:,:,1:end-1,:))/h;