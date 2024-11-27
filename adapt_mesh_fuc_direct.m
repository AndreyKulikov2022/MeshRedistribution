h=1;
N=40;
x=zeros(2,1,N,N);
[Q, K]=Q_field(N,N);
for i=1:N
    for j=1:N
        x(:,:,i,j)=[h*(j-1),h*(i-1)];
    end
end
dt=h.^2/2;
%metric=@(x) cartezian_metric(x);
%metric=@(x) rotated_metric(x);
metric=@(x) sin_metric(x,[1;0], N);
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
n_steps=40;
figure
hold on
h_q=plot([],[]);
B=zeros(size(x));

sub_set = zeros(1,1,size(x,3),size(x,4),"logical");
sub_set(:,:,1:2:end,1:2:end)=true;
sub_set(:,:,2:2:end,2:2:end)=true;

tau=1;
tol=0.2;

for i=1:n_steps
  x_new=gradient_decend(x,metric,Q, K);

  x_diff = x_new-x;

  %% Manually prevent singular mesh.
    tol=0.4;
    % Move half of the mesh and check if degeneration happened.
    x_temp = x + x_diff.*sub_set;
    % Check degeneration
    degen = is_degen(x_temp,tol);
    % Don't move to degenerate state.
    x_temp=x.*degen+x_temp.*(~degen);
    x=x_temp + x_diff.*(~sub_set);
    degen = is_degen(x,tol);
    x=x_temp.*degen+x.*(~degen);


    if video
        if draw_count==10
            delete(h_g);
            h_g=plot_mesh(x);
            pause(0.001);
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