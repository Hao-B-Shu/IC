%r=2, n=4    
       E=zeros(4); t=exp(1)^(2*pi*i/2);
       for k=1:4  E(k,k)=1; end
       A{1,1}=E(1,:)+E(2,:);   A{1,2}=E(3,:)+E(4,:);
       A{2,1}=E(1,:)+E(3,:);   A{2,2}=E(2,:)+E(4,:);
       B{1,1}=E(1,:)-E(2,:);   B{1,2}=E(3,:)-E(4,:);
       B{2,1}=E(1,:)-E(3,:);   B{2,2}=i*E(2,:)-i*E(4,:);
       C{1,1}=E(1,:)+i*E(2,:); C{1,2}=E(3,:)+i*E(4,:);
       C{2,1}=E(1,:)+i*E(3,:); C{2,2}=E(2,:)+i*E(4,:);
       D{1,1}=E(1,:)-i*E(2,:); D{1,2}=i*E(3,:)-i*i*E(4,:);
       D{2,1}=E(1,:)-i*E(3,:); D{2,2}=E(2,:)-i*E(4,:);
       for s=1:2
       for k=1:1  F{s,k}=B{s,1}+(t^(k-1))*B{s,2}; end
       end
       for s=1:2
       for k=1:1  G{s,k}=D{s,1}+(t^(k-1))*D{s,2}; end
       end
       l=1;  H=zeros(16);
       for s=1:2
       for k=1:2  H(l,:)=reshape(A{s,k}'*A{s,k},1,16); l=l+1; end
       end
       for s=1:2
       for k=1:2  H(l,:)=reshape(C{s,k}'*C{s,k},1,16); l=l+1; end
       end
       for s=1:2
       for k=1:1  H(l,:)=reshape(F{s,k}'*F{s,k},1,16); l=l+1; end
       end
       for s=1:2
       for k=1:1  H(l,:)=reshape(G{s,k}'*G{s,k},1,16); l=l+1; end
       end
       for s=1:4  H(l,:)=reshape(E(s,:)'*E(s,:),1,16); l=l+1; end
       rank(H)