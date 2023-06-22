% r>=5 might be changed.
r=5;
n=2*r;
t=exp(1)^(2*pi*i/r);
q=exp(1)^(2*pi*i/n);
E=zeros(n);
for k=1:n
E(k,k)=1;
end

for k=1:r
    for s=1:k
    A{2*k-1,s}=E(s,:)+E(2*k+1-s,:);
    A{2*k-1,s}=A{2*k-1,s}/norm(A{2*k-1,s});
    end
    if k<r-1
       for s=(k+1):(r-1)
       A{2*k-1,s}=E(s+k,:)+E(2*r+k-s,:);
       A{2*k-1,s}=A{2*k-1,s}/norm(A{2*k-1,s});
       end
    else
    end
    A{2*k-1,r}=E(r+k,:)+E(2*r,:);
    A{2*k-1,r}=A{2*k-1,r}/norm(A{2*k-1,r});
end

for k=1:r-1
    for s=1:k
    A{2*k,s}=E(s,:)+E(2*k-s+2,:);
    A{2*k,s}=A{2*k,s}/norm(A{2*k,s});
    end
    if k<r-1
       for s=(k+2):r
       A{2*k,s}=E(s+k,:)+E(2*r+k-s+1,:);
       A{2*k,s}=A{2*k,s}/norm(A{2*k,s});
       end
    else
    end
    A{2*k,k+1}=E(k+1,:)+E(2*r,:);
    A{2*k,k+1}=A{2*k,k+1}/norm(A{2*k,k+1});
end

for k=1:r
    for s=1:k
    B{2*k-1,s}=E(s,:)-E(2*k+1-s,:);
    end
    if k<r-1
       for s=(k+1):(r-1)
       B{2*k-1,s}=E(s+k,:)-E(2*r+k-s,:);
       end
    else
    end
    B{2*k-1,r}=E(r+k,:)-E(2*r,:);
end

for k=1:r-1
    for s=1:k
    B{2*k,s}=E(s,:)-E(2*k-s+2,:);
    end
    if k<r-1
       for s=(k+2):r
       B{2*k,s}=E(s+k,:)-E(2*r+k-s+1,:);
       end
    else
    end
    B{2*k,k+1}=E(k+1,:)-E(2*r,:);
end

for k=1:r
    for s=1:k
    C{2*k-1,s}=E(s,:)+q*E(2*k+1-s,:);
    C{2*k-1,s}=C{2*k-1,s}/norm(C{2*k-1,s});
    end
    if k<r-1
       for s=(k+1):(r-1)
       C{2*k-1,s}=E(s+k,:)+q*E(2*r+k-s,:);
       C{2*k-1,s}=C{2*k-1,s}/norm(C{2*k-1,s});
       end
    else
    end
    C{2*k-1,r}=E(r+k,:)+q*E(2*r,:);
    C{2*k-1,r}=C{2*k-1,r}/norm(C{2*k-1,r});
end

for k=1:r-1
    for s=1:k
    C{2*k,s}=E(s,:)+q*E(2*k-s+2,:);
    C{2*k,s}=C{2*k,s}/norm(C{2*k,s});
    end
    if k<r-1
       for s=(k+2):r
       C{2*k,s}=E(s+k,:)+q*E(2*r+k-s+1,:);
       C{2*k,s}=C{2*k,s}/norm(C{2*k,s});
       end
    else
    end
    C{2*k,k+1}=E(k+1,:)+q*E(2*r,:);
    C{2*k,k+1}=C{2*k,k+1}/norm(C{2*k,k+1});
end

for k=1:r
    for s=1:k
    D{2*k-1,s}=E(s,:)-q*E(2*k+1-s,:);
    end
    if k<r-1
       for s=(k+1):(r-1)
       D{2*k-1,s}=E(s+k,:)-q*E(2*r+k-s,:);
       end
    else
    end
    D{2*k-1,r}=E(r+k,:)-q*E(2*r,:);
end

for k=1:r-1
    for s=1:k
    D{2*k,s}=E(s,:)-q*E(2*k-s+2,:);
    end
    if k<r-1
       for s=(k+2):r
       D{2*k,s}=E(s+k,:)-q*E(2*r+k-s+1,:);
       end
    else
    end
    D{2*k,k+1}=E(k+1,:)-q*E(2*r,:);
end

F=cell(2*r-1,r);
for s=1:2*r-1
    for l=1:r
    F{s,l}=zeros(1,n);
    end
end

G=cell(2*r-1,r);
for s=1:2*r-1
    for l=1:r
    G{s,l}=zeros(1,n);
    end
end

for s=1:2*r-1
    for l=1:r
        for k=1:r
        F{s,l}=F{s,l}+(t^(k*l-k-l+1))*B{s,k};
        end
        F{s,l}=F{s,l}/norm(F{s,l});
    end
end

for s=1:2*r-1
    for l=1:r
        for k=1:r
        G{s,l}=G{s,l}+(t^(k*l-k-l+1))*D{s,k};
        end
        G{s,l}=G{s,l}/norm(G{s,l});
    end
end

for s=1:2*r-1
    for k=1:r
    A{s,k+r}=F{s,k};
    end
end

for s=1:2*r-1
    for k=1:r
    C{s,k+r}=G{s,k};
    end
end

for s=1:2*r-1
    for k=1:2*r
    A{s+2*r-1,k}=C{s,k};
    end
end

H=[];
for k=n+r-2:2*n-3
    Q=[];
    for s=1:n
    Q=[Q,reshape(A{k,s}'*A{k,s}-(1/n)*E,n*n,1)];
    end
    H=[H,orth(Q)];
end
Q=[];
for s=1:n
    Q=[Q,reshape(E(s,:)'*E(s,:)-(1/n)*E,n*n,1)];
end
H=[H,orth(Q)];
for g=1:r
    Q=[];
    for s=1:n
    Q=[Q,reshape(A{g,s}'*A{g,s}-(1/n)*E,n*n,1)];
    end
    H=[H,orth(Q)];
end
H=[H,zeros(n*n,1)];

d=eig(H);
f=d(2:end);
g=abs(prod(f,'all'))
