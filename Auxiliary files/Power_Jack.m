function [Jf]=Power_Jack(x)
global data
n=length(data.H);
D=data.D;
H=data.H;
%A=data.A;
gamma=data.gamma;
K=data.K;
omega=data.omega_R;
nn=length(x);
k=nn/2;
Jf = zeros(nn,nn); 
%f(1:nnodes)=x(nnodes+1:end); %x1dot =x2;
delta = x(1:k);   %x1
cal=zeros(n,n);
for i =1:n
  for j =1:n
     if i==j
         cal(i,j)=0;
     else
         cal(i,j)=K(i,j)*cos(delta(i)-delta(j)-gamma(i,j));
     end
  end
end
fi=sum(cal,2);
Jf(1:k,1:k)= - diag(D./(2*H)) - diag(omega./(2*H).*fi);

end