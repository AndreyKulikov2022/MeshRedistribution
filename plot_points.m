function h=plot_points(x,sub_set,color)
if(nargin<2)
    color = 'r';
end
sub_set=squeeze(sub_set);
Points=zeros(2,0);
for i=1:size(x,3)
    indexes=find(sub_set(:,i));
    Ps=x(:,:,indexes,i);
    Ps=squeeze(Ps);
    Points=[Points,Ps];
end
h=plot(Points(1,:),Points(2,:),'*');