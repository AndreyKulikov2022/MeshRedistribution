function [cd1,cd2]=diff_central(object,h)
cd1 = zeros(size(object)-[0,0,1,2]);
cd2 = zeros(size(object)-[0,0,2,1]);

cd1(:,:,:,:)=(object(:,:,1:end-1,3:end)-object(:,:,1:end-1,1:end-2)+object(:,:,2:end,3:end)-object(:,:,2:end,1:end-2))/(4*h);
cd2(:,:,:,:)=(object(:,:,3:end,1:end-1)-object(:,:,1:end-2,1:end-1)+object(:,:,3:end,2:end)-object(:,:,1:end-2,2:end))/(4*h);
