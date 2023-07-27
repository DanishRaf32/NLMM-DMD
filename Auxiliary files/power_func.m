function [f]=power_func(x)
global data
n=length(data.H);
D=data.D;
H=data.H;
A=data.A;
gamma=data.gamma;
K=data.K;
omega=data.omega_R;
nnodes=numel(x)/2;      
f = zeros(size(x)); 
f(1:nnodes)=x(nnodes+1:end); %x1dot =x2;
delta = x(1:nnodes);   %x1
cal=zeros(n,n);
for i =1:n
  for j =1:n
     if i==j
         cal(i,j)=0;
     else
         cal(i,j)=K(i,j)*sin(delta(i)-delta(j)-gamma(i,j));
     end
  end
end
fi=sum(cal,2);
f(nnodes+1:end)= - D./(2*H).*f(1:nnodes) + omega./(2*H).*A -...
                 omega./(2*H).*fi;

end