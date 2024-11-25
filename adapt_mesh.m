h=0.025;
N=1/h+1;
x=zeros(2,1,N,N);
I=eye_field(N,N);
for i=1:N
    for j=1:N
        x(:,:,i,j)=[h*(j-1),h*(i-1)];
    end
end
dt=h.^2/2;
%metric=@(x) cartezian_metric(x);
%metric=@(x) rotated_metric(x);
metric=@(x) sin_metric(x,[1;0]);
%metric=@(x) circ_metric(x,[1;0],1/2);
video=false;
if video
    figure("Position",[0,0,1000,1000])
    hold on
    testGIF='Mesh.gif';
    F=getframe(gcf);
    im=frame2im(F);
    [imind,cm] = rgb2ind(im,128);
    imwrite(imind,cm,testGIF,'gif','DelayTime',0.0001, 'Loopcount',inf);
    draw_count=10;
    plot([0:0.01:1],1/2+1/4*sin(2*pi*[0:0.01:1]),'r');
    h_g=plot_mesh(x);
end
n_steps=1000;
%x=fsolve(@(x)mmpde_rhs(x,metric,I,h),x);
%x=non_conservative(x,metric,I,h,dt,n_steps);
bs=zeros(1,n_steps);
bders=zeros(1,n_steps);
as=zeros(1,n_steps);
aders=zeros(1,n_steps);
rhss=zeros(1,n_steps);
figure
hold on
h_q=plot([],[]);
B=zeros(size(x));

sub_set = zeros(1,1,size(x,3),size(x,4),"logical");
sub_set(:,:,1:2:end,1:2:end)=true;
sub_set(:,:,2:2:end,2:2:end)=true;

for i=1:n_steps
  %  x=semi_conservative(x,metric,I,h,dt,n_steps);
    [rhs,bs(i),bders(i),as(i),aders(i),q_ali]=mmpde_rhs(x,metric,I,h);
    rhss(i)=sum(abs(rhs),"all");

    %% Manually prevent singular mesh.
    tol=0.4;
    % Move half of the mesh and check if degeneration happened.
    x_temp = x + dt*rhs.*sub_set;
    % Check degeneration
    degen = is_degen(x_temp,tol);
    % Don't move to degenerate state.
    x_temp=x.*degen+x_temp.*(~degen);
    x=x_temp + dt*rhs.*(~sub_set);
    degen = is_degen(x,tol);
    x=x_temp.*degen+x.*(~degen);

%     if(any(degen,"all"))
%         disp("oops");
%     end

   % x=x+dt*rhs;
%     delete(h_q);
%     h_q=plot4d(x,q_ali);
%     clim([0,30]);
%     colorbar
%     pause(0.01);
    if video
        if draw_count==10
            delete(h_g);
            h_g=plot_mesh(x);
            %pause(0.001);
            F=getframe(gcf);
            im=frame2im(F);
            [imind,cm] = rgb2ind(im,256);
            imwrite(imind,cm,testGIF,'gif','DelayTime',0.0001,'WriteMode','append');
            draw_count=0;
        else
            draw_count=draw_count+1;
        end
    end
end
figure("Position",[0,0,1000,1000])
hold on
plot([0:0.01:1],1/2+1/4*sin(2*pi*[0:0.01:1]),'r');
h_g=plot_mesh(x);
% figure
% hold on
% plot([1:n_steps],bs);
% title("B")
% figure
% hold on
% plot([1:n_steps],bders);
% title("B*derivative")
% figure
% hold on
% plot([1:n_steps],as);
% title("A")
% figure
% hold on
% plot([1:n_steps],aders);
% title("A*derivative")
% 
% figure
% hold on
% plot([1:n_steps],rhss);
% title("rhs")
% lb=zeros(size(x));
% ub=ones(size(x));
% ub(1,:,:,1)=lb(1,:,:,1);
% lb(1,:,:,end)=ub(1,:,:,end);
% ub(2,:,1,:)=lb(2,:,1,:);
% lb(2,:,end,:)=rb(2,:,end,:);
% %constrain jacobian
% %lsqnonlin(@myfun,x0,lb,ub,A,b,Aeq,beq,@mycon,options),@(x)Jconstr(x,h)
% A = [];
% b = [];
% Aeq = [];
% beq = [];
% x2=lsqnonlin(@(x)calc_q_ali(x,metric,h),x,lb,ub,A,b,Aeq,beq,@(x)Jconstr(x,h));