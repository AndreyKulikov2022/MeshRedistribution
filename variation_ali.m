function [result, p,bs,bders,as,aders,q_ali] = variation_ali(x,metric,I, h)
%% M, I - 4-D arrays, first two are for the matrix itself, third and second is for 2D domain 
% x - coodrdinates
% M - metric tensor, I - identity matrices.
% h - uniform grid step
result=zeros(size(x));
[d_xi1,d_xi2] = diff_xi(x,h);
J=cat(2,d_xi1,d_xi2);
Ji=pageinv(J);
d11 = zeros(size(x));
d12 = zeros(size(x));
d22 = zeros(size(x));
d11(:,:,:,2:end-1)=(x(:,:,:,3:end)-2*x(:,:,:,2:end-1)+x(:,:,:,1:end-2))/(h^2);
d22(:,:,2:end-1,:)=(x(:,:,3:end,:)-2*x(:,:,2:end-1,:)+x(:,:,1:end-2,:))/(h^2);
d12(:,:,2:end-1,2:end-1)=(x(:,:,3:end,3:end)-x(:,:,3:end,1:end-2)-x(:,:,1:end-2,3:end)+x(:,:,1:end-2,1:end-2))/(h^2);
% [d2_xi11,d2_xi12] = diff_xi(d_xi1,h);
% [d2_xi21,d2_xi22] = diff_xi(d_xi2,h);
%d2s={d2_xi11,d2_xi12;d2_xi21,d2_xi22};
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
% Balancing function
frob2=@(F) sum(sum(F.*F,1),2);
p=sqrt(frob2(As{1,1})+frob2(As{1,2})+frob2(As{2,1})+frob2(As{2,2})+frob2(Bs{1})+frob2(Bs{2}));
%% Euler-Lagrange and balancing function
%Convective term is high, depending on Bs sign different stencil for
%derivative might be needed
Js=cell(2,1);
d_xi1_central = zeros(size(x));
d_xi1_central(:,:,2:end-1,2:end-1)=(x(:,:,2:end-1,3:end)-x(:,:,2:end-1,1:end-2))/(2*h);
d_xi1_fb = zeros(size(x));
d_xi1_forward=(x(:,:,2:end-1,3:end)-x(:,:,2:end-1,2:end-1))/(h);
d_xi1_backward=(x(:,:,2:end-1,2:end-1)-x(:,:,2:end-1,1:end-2))/(h);
ids1=Bs{1}(1,1,:,:)>0;
d_xi1_fb(:,:,2:end-1,2:end-1)=d_xi1_forward.*ids1(:,:,2:end-1,2:end-1)+d_xi1_backward.*(~ids1(:,:,2:end-1,2:end-1));
w=frob2(Bs{1})./p;
Js{1}=d_xi1_fb;%.*w + d_xi1_central.*(1-w);
d_xi2_central = zeros(size(x));
d_xi2_central(:,:,2:end-1,2:end-1)=(x(:,:,3:end,2:end-1)-x(:,:,1:end-2, 2:end-1))/(2*h);
d_xi2_fb = zeros(size(x));
d_xi2_forward=(x(:,:,3:end,2:end-1)-x(:,:,2:end-1,2:end-1))/(h);
d_xi2_backward=(x(:,:,2:end-1,2:end-1)-x(:,:,1:end-2,2:end-1))/(h);
ids2=Bs{2}(1,1,:,:)>0;
d_xi2_fb(:,:,2:end-1,2:end-1)=d_xi2_forward.*ids2(:,:,2:end-1,2:end-1)+d_xi2_backward.*(~ids2(:,:,2:end-1,2:end-1));
w=frob2(Bs{2})./p;
Js{2}=d_xi2_fb;%.*w + d_xi2_central.*(1-w);
% More than that where difusive term is dominant central difference is
% better.
% coef_sum=abs(Bs{1}(1,1,:,:))+abs(Bs{2}(1,1,:,:))+abs(As{1,1}(1,1,:,:))/h+abs(As{2,2}(1,1,:,:))/h+abs(As{1,2}(1,1,:,:))/h+abs(As{2,1}(1,1,:,:))/h;
% w=abs(Bs{1}(1,1,:,:))./coef_sum;
% Js{1}=Js{1}.*w+d_xi1.*(1-w);
% w=abs(Bs{2}(1,1,:,:))./coef_sum;
% Js{2}=Js{2}.*w+d_xi2.*(1-w);
d2s={d11,d12;d12,d22};
bs=0;
bders=0;
as=0;
aders=0;
for i=1:2
    %result=result+pagemtimes(Bs{i},J(:,i,:,:));
    bs=max([bs,sum(abs(Bs{i}),"all")]);
    bders=max([bders,sum(abs(pagemtimes(Bs{i},Js{i})),"all")]);
    result=result+pagemtimes(Bs{i},Js{i});
    for j=1:2
        as=max([as,sum(abs(As{i,j}),"all")]);
        aders=max([aders,sum(abs(pagemtimes(As{i,j},d2s{i,j})),"all")]);
        result=result+pagemtimes(As{i,j},d2s{i,j});
    end
end
result(:,:,1,:)=0;
result(:,:,end,:)=0;
result(:,:,:,1)=0;
result(:,:,:,end)=0;

q_ali=zeros(1,1,size(x,3),size(x,4));
JMJ=pagemtimes(J,'transpose', ...
            pagemtimes(M,J),'none');
%makes sense, make the stencil variable based on Bs and As, Bs are bigger,
%but nor everywhere.
every=abs(Bs{1}(1,1,:,:))+abs(Bs{2}(1,1,:,:))+abs(As{1,1}(1,1,:,:))+abs(As{2,2}(1,1,:,:));
% every=abs(pagemtimes(Bs{1},Js{1}))+abs(pagemtimes(Bs{2},Js{2}))+abs(pagemtimes(As{1,1},d2s{1,1}))+abs(pagemtimes(As{2,2},d2s{2,2}));
% temp=pagemtimes(Bs{1},Js{1});
q_ali=(JMJ(1,1,:,:)+JMJ(2,2,:,:))./(2*sqrt(JMJ(1,1,:,:).*JMJ(2,2,:,:)-JMJ(1,2,:,:).*JMJ(2,1,:,:)));
end