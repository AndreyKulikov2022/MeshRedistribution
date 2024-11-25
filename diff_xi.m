function [d_xi1,d_xi2] = diff_xi(object,h)
d_xi1 = zeros(size(object));
d_xi2 = zeros(size(object));

%d_xi1(:,:,:,2:end-1)=(object(:,:,:,3:end)-object(:,:,:,2:end-1))/(h);
d_xi1(:,:,:,2:end-1)=(object(:,:,:,3:end)-object(:,:,:,1:end-2))/(2*h);
d_xi1(:,:,:,1)=(-3*object(:,:,:,1)+4*object(:,:,:,2) -object(:,:,:,3))/(2*h);
d_xi1(:,:,:,end)=-(-3*object(:,:,:,end)+4*object(:,:,:,end-1) -object(:,:,:,end-2))/(2*h);

%d_xi2(:,:,2:end-1,:)=(object(:,:,3:end,:)-object(:,:,2:end-1,:))/(h);
d_xi2(:,:,2:end-1,:)=(object(:,:,3:end,:)-object(:,:,1:end-2,:))/(2*h);
d_xi2(:,:,1,:)=(-3*object(:,:,1,:)+4*object(:,:,2,:) -object(:,:,3,:))/(2*h);
d_xi2(:,:,end,:)=-(-3*object(:,:,end,:)+4*object(:,:,end-1,:) -object(:,:,end-2,:))/(2*h);

end

% Do honest implicit Euler