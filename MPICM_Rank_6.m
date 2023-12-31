%r=3, n=6	
    E=zeros(6);  t=exp(1)^(2*pi*i/3);
	for k=1:6  E(k,k)=1;  end
	A{1,1}=E(1,:)+E(2,:);  A{1,2}=E(3,:)+E(5,:);  A{1,3}=E(4,:)+E(6,:);
	A{2,1}=E(1,:)+E(3,:);  A{2,2}=E(2,:)+E(6,:);  A{2,3}=E(4,:)+E(5,:);
	A{3,1}=E(1,:)+E(4,:);  A{3,2}=E(2,:)+E(3,:);  A{3,3}=E(5,:)+E(6,:);
	B{1,1}=E(1,:)-E(2,:);   B{1,2}=E(3,:)-E(5,:);   B{1,3}=i*E(4,:)-i*E(6,:);
	B{2,1}=E(1,:)-E(3,:);   B{2,2}=E(2,:)-E(6,:);   B{2,3}=E(4,:)-E(5,:);
	B{3,1}=E(1,:)-E(4,:);   B{3,2}=E(2,:)-E(3,:);   B{3,3}=E(5,:)-E(6,:);
	C{1,1}=E(1,:)+i*E(2,:);  C{1,2}=E(3,:)+i*E(5,:);  C{1,3}=E(4,:)+i*E(6,:);
	C{2,1}=E(1,:)+i*E(3,:);  C{2,2}=E(2,:)+i*E(6,:);  C{2,3}=E(4,:)+i*E(5,:);
	C{3,1}=E(1,:)+i*E(4,:);  C{3,2}=E(2,:)+i*E(3,:);  C{3,3}=E(5,:)+i*E(6,:);
	D{1,1}=E(1,:)-i*E(2,:);   D{1,2}=E(3,:)-i*E(5,:);   D{1,3}=E(4,:)-i*E(6,:);
	D{2,1}=E(1,:)-i*E(3,:);   D{2,2}=E(2,:)-i*E(6,:);   D{2,3}=E(4,:)-i*E(5,:);
	D{3,1}=E(1,:)-i*E(4,:);   D{3,2}=E(2,:)-i*E(3,:);   D{3,3}=E(5,:)-i*E(6,:);
	for s=1:3
	for k=1:2  F{s,k}=B{s,1}+(t^(k-1))*B{s,2}+(t^(2*k-2))*B{s,3}; end
	end
	for s=1:3
	for k=1:2  G{s,k}=D{s,1}+(t^(k-1))*D{s,2}+(t^(2*k-2))*D{s,3}; end
	end
	l=1;  H=zeros(36);
	for s=1:3
	for k=1:3  H(l,:)=reshape(A{s,k}'*A{s,k},1,36);  l=l+1; end
	end
	for s=1:3
	for k=1:3  H(l,:)=reshape(C{s,k}'*C{s,k},1,36);  l=l+1; end
	end
	for s=1:3
	for k=1:2  H(l,:)=reshape(F{s,k}'*F{s,k},1,36);  l=l+1;  end
	end
	for s=1:3
	for k=1:2  H(l,:)=reshape(G{s,k}'*G{s,k},1,36);  l=l+1; end
	end
	for s=1:6      H(l,:)=reshape(E(s,:)'*E(s,:),1,36);     l=l+1;  end
	rank(H)