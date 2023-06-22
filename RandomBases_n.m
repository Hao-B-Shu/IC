% n might be changed
n=10;
R=cell(1,n+1);
I=eye(n);

for k=1:n+1
    r1=rand(1,n^2);
    r2=rand(1,n^2);
    r=r1+r2*i;
    R{1,k}=orth(reshape(r,n,n));
    if rank(R{1,k})~=n
        rank(R{1,k})
    end
end

for l=1:n+1
    Q=[];
    for k=1:n
        Q=[Q,reshape(R{1,l}(:,k)*R{1,l}(:,k)'/(norm(R{1,l}(:,k))^2)-(1/n)*I,n*n,1)];
    end
    R{1,l}=orth(Q);
end

H=[];
for l=1:n+1
    H=[H,R{1,l}];
end
H=[H,zeros(n*n,1)];

d=eig(H);
f=d(2:end);
g=abs(prod(f,'all'))
