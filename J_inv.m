function Ji = J_inv(I, x,h)
%% There are 4 zones
    %       (j,k+1)
    %          |
    %        1 | 2
    %(j-1,k)___|___(j+1,k)
    %        3 | 4
    %          |
    %       (j,k-1)
% Must find derivatives in each and average.
x_jk = x(:,:,2:end-1,2:end-1);
x_jp_k = x(:,:,2:end-1,3:end);
x_jm_k = x(:,:,2:end-1,1:end-2);
x_j_kp = x(:,:,3:end,2:end-1);
x_j_km = x(:,:,1:end-2,2:end-1);
Ji=(one_zone(I,h, x_jm_k,x_jk,x_jk,x_j_kp) + one_zone(I,h, x_jk,x_jp_k,x_jk,x_j_kp) +...
    one_zone(I,h, x_jm_k,x_jk,x_j_km,x_jk) + one_zone(I,h, x_jk,x_jp_k,x_j_km,x_jk))/4;
% J(:,1,[1,size(x,3)],2:end-1)=(x(:,:,[1,size(x,3)],3:end)-x(:,:,[1,size(x,3)],1:end-2))/(2*h);
% J(:,2,1,:)=(x(:,:,2,:)-x(:,:,1,:))/h;
% J(:,2,end,:)=(x(:,:,end,:)-x(:,:,end-1,:))/h;
% 
% J(:,2,2:end-1,[1,size(x,4)])=(x(:,:,3:end,[1,size(x,4)])-x(:,:,1:end-2,[1,size(x,4)]))/(2*h);
% J(:,1,:,1)=(x(:,:,:,2)-x(:,:,:,1))/h;
% J(:,1,:,end)=(x(:,:,:,end)-x(:,:,:,end-1))/h;
end

function Ji = one_zone(I,h, x1,x2,x3,x4)
J=I;
J(:,1,2:end-1,2:end-1)=(x2-x1)/h;
J(:,2,2:end-1,2:end-1)=(x4-x3)/h;
Ji=pageinv(J);
end