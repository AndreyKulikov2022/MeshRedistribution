function x_new = implicit_euler(x,metric,I, h, tau, dt)
%% Right hand side coefficients in the Euler-Lagrange equations
% Boundary components are zeroed out.
[As,Bs] = A_B_ali(x,metric,I, h);
%% Balance
frob2=@(F) sum(sum(F.*F,1),2);
p=sqrt(frob2(As{1,1})+frob2(As{1,2})+frob2(As{2,1})+frob2(As{2,2})+frob2(Bs{1})+frob2(Bs{2}));
As{1,1}=As{1,1}./p/tau/h^2*dt;
As{1,2}=As{1,2}./p/tau/4/h^2*dt;
As{2,1}=As{2,1}./p/tau/4/h^2*dt;
As{2,2}=As{2,2}./p/tau/h^2*dt;
Bs{1}=Bs{1}./p/tau/(h)*dt;
Bs{2}=Bs{2}./p/tau/(h)*dt;

 %% Zero boundary components
           As{1,1}(:,:,1,:) = 0;
           As{1,1}(:,:,end,:) = 0;
           As{1,1}(:,:,:,1) = 0;
           As{1,1}(:,:,:,end) = 0;

           As{1,2}(:,:,1,:) = 0;
           As{1,2}(:,:,end,:) = 0;
           As{1,2}(:,:,:,1) = 0;
           As{1,2}(:,:,:,end) = 0;

           As{2,1}(:,:,1,:) = 0;
           As{2,1}(:,:,end,:) = 0;
           As{2,1}(:,:,:,1) = 0;
           As{2,1}(:,:,:,end) = 0;

           As{2,2}(:,:,1,:) = 0;
           As{2,2}(:,:,end,:) = 0;
           As{2,2}(:,:,:,1) = 0;
           As{2,2}(:,:,:,end) = 0;

           Bs{1}(:,:,1,:) = 0;
           Bs{1}(:,:,end,:) = 0;
           Bs{1}(:,:,:,1) = 0;
           Bs{1}(:,:,:,end) = 0;

           Bs{2}(:,:,1,:) = 0;
           Bs{2}(:,:,end,:) = 0;
           Bs{2}(:,:,:,1) = 0;
           Bs{2}(:,:,:,end) = 0;   
        %

%% Squeeze x, y, As and Bs into 1D vectors and form the matrix.
% tau * p/ dt * [x_{n+1}-x_n; y_{n+1}-y_n] = A * [x^{n+1}; y^{n+1}];
% [x_{n+1}-x_n; y_{n+1}-y_n] = (A11*dxi2 + (A12+A21)dxideta + A22deta2 + B1dxi + B2deta)*dt/(tau * p) * [x^{n+1}; y^{n+1}];
% (eye(N*N) - (A11*dxi2 + (A12+A21)dxideta + A22deta2 + B1dxi + B2deta)*dt/(tau * p))[x_{n+1}; y_{n+1}] =[x^{n}; y^{n}];
%Explicitely use the fact that in this example A and B are diagonal.
%Otherwise one more matrix multiplication will be needed.
%squeezing will be done by pairs (xi,yi), xs are odd, ys are even
N=size(x,3);
x_1d = zeros(2*N^2,1);
A11 = zeros(2*N^2,1);
A12 = zeros(2*N^2,1);
A22 = zeros(2*N^2,1);
B1 = zeros(2*N^2,1);
B2 = zeros(2*N^2,1);
for i=0:N-1
x_1d(2*i*N+1:2:2*(i+1)*N)=x(1,1,i+1,:);
x_1d(2*i*N+2:2:2*(i+1)*N)=x(2,1,i+1,:);

A11(2*i*N+1:2:2*(i+1)*N)=As{1,1}(1,1,i+1,:);
A11(2*i*N+2:2:2*(i+1)*N)=As{1,1}(2,2,i+1,:);

A12(2*i*N+1:2:2*(i+1)*N)=As{1,2}(1,1,i+1,:)+As{2,1}(1,1,i+1,:);
A12(2*i*N+2:2:2*(i+1)*N)=As{1,2}(2,2,i+1,:)+As{2,1}(2,2,i+1,:);

A22(2*i*N+1:2:2*(i+1)*N)=As{2,2}(1,1,i+1,:);
A22(2*i*N+2:2:2*(i+1)*N)=As{2,2}(2,2,i+1,:);

B1(2*i*N+1:2:2*(i+1)*N)=Bs{1}(1,1,i+1,:);
B1(2*i*N+2:2:2*(i+1)*N)=Bs{1}(2,2,i+1,:);

B2(2*i*N+1:2:2*(i+1)*N)=Bs{2}(1,1,i+1,:);
B2(2*i*N+2:2:2*(i+1)*N)=Bs{2}(2,2,i+1,:);
end
   
A = eye(2*(N^2));
ids1=B1>0;%forward
ids2=B2>0;
%% Filling the matrix
% Point exactly - at the main diagonal,
% Neighbours in xi direction - on the diag(2)
% Neighbours in eta seem to be at diag(2*N)
% A = A-(diag(-2*(A11+A22))+diag(A11(1:end-2)+B1(1:end-2),2)+diag(A11(3:end)-B1(3:end),-2)...
%     +diag(A22(1:end-2*N)+B2(1:end-2*N),2*N)+diag(A22(2*N+1:end)-B2(2*N+1:end),-2*N)...
%     +diag(A12(1:end-2*N-2),2*N+2)+diag(-A12(1:end-2*(N-1)),2*N-2)+diag(A12(2*N+3:end),-2*N-2)+diag(-A12(2*N-1:end),-2*N+2));
A = A-(diag(-2*(A11+A22))+diag(A11(1:end-2),2)+diag(A11(3:end),-2)...
    +diag(A22(1:end-2*N),2*N)+diag(A22(2*N+1:end),-2*N)...
    +diag(A12(1:end-2*N-2),2*N+2)+diag(-A12(1:end-2*(N-1)),2*N-2)+diag(A12(2*N+3:end),-2*N-2)+diag(-A12(2*N-1:end),-2*N+2)...
    +diag(-B1(3:end).*(~ids1(3:end)),-2)+diag(-B1.*ids1 + B1.*(~ids1))+diag(B1(1:end-2).*ids1(1:end-2),2)...
    +diag(-B2(2*N+1:end).*(~ids2(2*N+1:end)),-2*N)+diag(-B2.*ids2 + B2.*(~ids2))+diag(B2(1:end-2*N).*ids2(1:end-2*N),2*N));


x_1d_new=linsolve(A,x_1d);
x_new=x;
for i=0:N-1
x_new(1,1,i+1,:)=x_1d_new(2*i*N+1:2:2*(i+1)*N);
x_new(2,1,i+1,:)=x_1d_new(2*i*N+2:2:2*(i+1)*N);
end
