function h=plot_mesh(x, color)
if(nargin<2)
    color = 'b';
end
Points=zeros(2,0);
for i=1:size(x,3)
    Ps=x(:,:,:,i);
    Ps=squeeze(Ps);
    Points=[Points,Ps,[nan;nan]];
    Ps=x(:,:,i,:);
    Ps=squeeze(Ps);
    Points=[Points,Ps,[nan;nan]];
end
h=plot(Points(1,:),Points(2,:),color);