function [d_xi1_fb,d_xi2_fb] = diff_xi_fb(object,h,B)
x=object;

d_xi1_fb = zeros(size(x));
d_xi1_forward=(x(:,:,2:end-1,3:end)-x(:,:,2:end-1,2:end-1))/(h);
d_xi1_backward=(x(:,:,2:end-1,2:end-1)-x(:,:,2:end-1,1:end-2))/(h);
ids1=B(1,1,:,:)>0;
d_xi1_fb(:,:,2:end-1,2:end-1)=d_xi1_forward.*ids1(:,:,2:end-1,2:end-1)+d_xi1_backward.*(~ids1(:,:,2:end-1,2:end-1));

d_xi2_fb = zeros(size(x));
d_xi2_forward=(x(:,:,3:end,2:end-1)-x(:,:,2:end-1,2:end-1))/(h);
d_xi2_backward=(x(:,:,2:end-1,2:end-1)-x(:,:,1:end-2,2:end-1))/(h);
ids2=B(2,1,:,:)>0;
d_xi2_fb(:,:,2:end-1,2:end-1)=d_xi2_forward.*ids2(:,:,2:end-1,2:end-1)+d_xi2_backward.*(~ids2(:,:,2:end-1,2:end-1));

d_xi1_fb(:,:,:,1)=(-3*object(:,:,:,1)+4*object(:,:,:,2) -object(:,:,:,3))/(2*h);
d_xi1_fb(:,:,:,end)=-(-3*object(:,:,:,end)+4*object(:,:,:,end-1) -object(:,:,:,end-2))/(2*h);

d_xi2_fb(:,:,1,:)=(-3*object(:,:,1,:)+4*object(:,:,2,:) -object(:,:,3,:))/(2*h);
d_xi2_fb(:,:,end,:)=-(-3*object(:,:,end,:)+4*object(:,:,end-1,:) -object(:,:,end-2,:))/(2*h);

end


