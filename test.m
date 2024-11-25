%unit square
h=0.1;
N=1/h+1;
x=zeros(2,1,N,N);
I=eye_field(N,N);
for i=1:N
    for j=1:N
        x(:,:,i,j)=[h*(j-1),h*(i-1)];
    end
end
[d_xi1,d_xi2] = diff_xi(x,h);
J=cat(2,d_xi1,d_xi2);
assert(all(all(all(all(round(J,6)==I)))));

%Evaluate metric
metric=@(x) cartezian_metric(x);
M=zeros(2,2,size(x,3),size(x,4));
ro=zeros(1,1,size(x,3),size(x,4));
for i=1:size(x,3)
    for j=1:size(x,4)
        M(:,:,i,j)=metric(x(:,:,i,j));
        ro(:,:,i,j)=det(M(:,:,i,j));
    end
end
assert(all(all(all(all(M==I)))));
[result, p] = variation_ali(x,metric,I, h);