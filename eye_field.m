function I=eye_field(n,m)
I=zeros(2,2,n,m);
for i=1:n
    for j=1:m
        I(:,:,i,j)=eye(2,2);
    end
end
end