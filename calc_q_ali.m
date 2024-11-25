function q_ali=calc_q_ali(x,metric, h)
[d_xi1,d_xi2] = diff_xi(x,h);
J=cat(2,d_xi1,d_xi2);
M=zeros(2,2,size(x,3),size(x,4));
for i=1:size(x,3)
    for j=1:size(x,4)
        M(:,:,i,j)=metric(x(:,:,i,j));
    end
end
JMJ=pagemtimes(J,'transpose', ...
            pagemtimes(M,J),'none');
q_ali=(JMJ(1,1,:,:)+JMJ(2,2,:,:))./(2*sqrt(JMJ(1,1,:,:).*JMJ(2,2,:,:)-JMJ(1,2,:,:).*JMJ(2,1,:,:)));
end