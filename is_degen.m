function degen = is_degen(x,tol)
degen = zeros(1,1,size(x,3),size(x,4),"logical");
v1=x(:,:,2:end-1,1:end-2) - x(:,:,2:end-1,2:end-1);
v2=x(:,:,1:end-2,2:end-1) - x(:,:,2:end-1,2:end-1);
v3=x(:,:,2:end-1,3:end) - x(:,:,2:end-1,2:end-1);
v4=x(:,:,3:end,2:end-1) - x(:,:,2:end-1,2:end-1);

v1=v1./sqrt(sum(v1.*v1,1));
v2=v2./sqrt(sum(v2.*v2,1));
v3=v3./sqrt(sum(v3.*v3,1));
v4=v4./sqrt(sum(v4.*v4,1));

%Find crosss products to monitor angles.
c1=v1(1,1,:,:).*v2(2,1,:,:)-v1(2,1,:,:).*v2(1,1,:,:);
c2=v2(1,1,:,:).*v3(2,1,:,:)-v2(2,1,:,:).*v3(1,1,:,:);
c3=v3(1,1,:,:).*v4(2,1,:,:)-v3(2,1,:,:).*v4(1,1,:,:);
c4=v4(1,1,:,:).*v1(2,1,:,:)-v4(2,1,:,:).*v1(1,1,:,:);

degen(:,:,2:end-1,2:end-1)=(c1<tol)|(c2<tol)|(c3<tol)|(c4<tol);
degen(:,:,2:end-1,2:end-1)=degen(:,:,2:end-1,2:end-1)|degen(:,:,1:end-2,2:end-1)|...
    degen(:,:,3:end,2:end-1)|degen(:,:,2:end-1,1:end-2)|degen(:,:,2:end-1,3:end);