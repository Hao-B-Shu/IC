%For n=2r, r in 5:50, calculate the ranks of operators derived by an MPICM
%Can calculate for a single r>=5 by replacing the first loop by a special r 
for r=5:50
	n=2*r;  t=exp(1)^(2*pi*i/r);  q=exp(1)^(2*pi*i/n);  E=zeros(n);
	for k=1:n    E(k,k)=1;  end
	for k=1:r
	for s=1:k  A{2*k-1,s}=E(s,:)+E(2*k+1-s,:);  end
	if k<r-1
	for s=(k+1):(r-1)  A{2*k-1,s}=E(s+k,:)+E(2*r+k-s,:);  end
	else
	end
	A{2*k-1,r}=E(r+k,:)+E(2*r,:);
	end
	for k=1:r-1
	for s=1:k  A{2*k,s}=E(s,:)+E(2*k-s+2,:); end
	if k<r-1
	for s=(k+2):r  A{2*k,s}=E(s+k,:)+E(2*r+k-s+1,:);  end
	else
	end
	A{2*k,k+1}=E(k+1,:)+E(2*r,:);
	end
	for k=1:r
	for s=1:k  B{2*k-1,s}=E(s,:)-E(2*k+1-s,:); end
	if k<r-1
	for s=(k+1):(r-1)  B{2*k-1,s}=E(s+k,:)-E(2*r+k-s,:); end
	else
	end
	B{2*k-1,r}=E(r+k,:)-E(2*r,:);
	end
	for k=1:r-1
	for s=1:k  B{2*k,s}=E(s,:)-E(2*k-s+2,:); end
	if k<r-1
	for s=(k+2):r  B{2*k,s}=E(s+k,:)-E(2*r+k-s+1,:); end
	else
	end
	B{2*k,k+1}=E(k+1,:)-E(2*r,:);
	end
	for k=1:r
	for s=1:k  C{2*k-1,s}=E(s,:)+q*E(2*k+1-s,:); end
	if k<r-1
	for s=(k+1):(r-1)  C{2*k-1,s}=E(s+k,:)+q*E(2*r+k-s,:); end
	else
	end
	C{2*k-1,r}=E(r+k,:)+q*E(2*r,:);
	end
	for k=1:r-1
	for s=1:k  C{2*k,s}=E(s,:)+q*E(2*k-s+2,:); end
	if k<r-1
	for s=(k+2):r  C{2*k,s}=E(s+k,:)+q*E(2*r+k-s+1,:); end
	else
	end
	C{2*k,k+1}=E(k+1,:)+q*E(2*r,:);
	end
	for k=1:r
	for s=1:k  D{2*k-1,s}=E(s,:)-q*E(2*k+1-s,:); end
	if k<r-1
	for s=(k+1):(r-1)  D{2*k-1,s}=E(s+k,:)-q*E(2*r+k-s,:); end
	else
	end
	D{2*k-1,r}=E(r+k,:)-q*E(2*r,:);
	end
	for k=1:r-1
	for s=1:k  D{2*k,s}=E(s,:)-q*E(2*k-s+2,:); end
	if k<r-1
	for s=(k+2):r  D{2*k,s}=E(s+k,:)-q*E(2*r+k-s+1,:); end
	else
	end
	D{2*k,k+1}=E(k+1,:)-q*E(2*r,:);
	end
	F=cell(2*r-1,r-1);
	for s=1:2*r-1
	for l=1:(r-1)  F{s,l}=zeros(1,n); end
	end
	G=cell(2*r-1,r-1);
	for s=1:2*r-1
	for l=1:(r-1)  G{s,l}=zeros(1,n); end
	end
	for s=1:2*r-1
	for l=1:(r-1)
	for k=1:r  F{s,l}=F{s,l}+(t^(k*l-k-l+1))*B{s,k}; end
	end
	end
	for s=1:2*r-1
	for l=1:(r-1)
	for k=1:r  G{s,l}=G{s,l}+(t^(k*l-k-l+1))*D{s,k}; end
	end
	end
	for s=1:2*r-1
	for k=1:r-1  A{s,k+r}=F{s,k}; end
	end
	for s=1:2*r-1
	for k=1:r-1  C{s,k+r}=G{s,k}; end
	end
	for s=1:2*r-1
	for k=1:2*r-1  A{s+2*r-1,k}=C{s,k}; end
	end
	H=zeros(n*n);  l=1;
	for k=n+r-2:2*n-3
	for s=1:n-1  H(l,:)=reshape(A{k,s}'*A{k,s},1,n*n);  l=l+1; end
	end
	for s=1:n  H(l,:)=reshape(E(s,:)'*E(s,:),1,n*n);  l=l+1; end
	for g=1:r
	for s=1:n-1  H(l,:)=reshape(A{g,s}'*A{g,s},1,n*n);  l=l+1; end
	end
	disp(['n=2r=',num2str(n)])
	disp(['Rank=',num2str(rank(H))])
end
