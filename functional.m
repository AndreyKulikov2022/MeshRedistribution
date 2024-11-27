function F = functional(x, metric, Q, K)

M=zeros(2,2,size(x,3),size(x,4));
ro=zeros(1,1,size(x,3),size(x,4));

for i=1:size(x,3)
    for j=1:size(x,4)
        M(:,:,i,j)=metric(x(:,:,i,j));
        ro(:,:,i,j)=sqrt(det(M(:,:,i,j)));
    end
end

%% Cell
%
% .(X01, Y01) .(X11, Y11)
%
% .(X00, YOO) .(X10, Y1O)

X00=x(1,1,1:end-1,1:end-1);
X10=x(1,1,1:end-1,2:end);
X01=x(1,1,2:end, 1:end-1);
X11=x(1,1,2:end, 2:end);
Xc = cat(1,X00,X10,X01,X11);

Y00=x(2,1,1:end-1,1:end-1);
Y10=x(2,1,1:end-1,2:end);
Y01=x(2,1,2:end, 1:end-1);
Y11=x(2,1,2:end, 2:end);
Yc = cat(1,Y00,Y10,Y01,Y11);

Ms={M(:,:,1:end-1,1:end-1),M(:,:,1:end-1,2:end),M(:,:,2:end, 1:end-1),M(:,:,2:end, 2:end)};
ros={ro(:,:,1:end-1,1:end-1),ro(:,:,1:end-1,2:end),ro(:,:,2:end, 1:end-1),ro(:,:,2:end, 2:end)};

F=zeros(size(X00));
for i=1:4
a1 = pagemtimes(Q{i},Xc);
a2 = pagemtimes(Q{i},Yc);
J = pagemtimes(a1,'transpose', ...
            pagemtimes(K,a2),'none');
F=F+(Ms{i}(1,1,:,:).*pagemtimes(a1,'transpose',a1,'none')...
    +2*Ms{i}(1,2,:,:).*pagemtimes(a1,'transpose',a2,'none')...
    +Ms{i}(2,2,:,:).*pagemtimes(a2,'transpose',a2,'none'))./J./ros{i};
end
F=sum(F,"all");