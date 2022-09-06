function [Opt_rank]=Consensus(X,Type,Multi_type)
%%%%%%%%%%%%%%%%%%%%
%A function for consensus-based community detection using both single-
%and multi-layer Louvain algorthim.

%Input:
        %X: The multiplex network
        %Community detection Type: Monolayer=1
                            %      Multilayer=2
        %Multi_type: No. of Participants/Tasks to segment multiplex network
        
        
%Output:
        %Opt_rank: The optimal number of clusters to extract
%%%%%%%%%%%%%%%%%%%%%

if Type==1
    Multi_type=0;
    Vs=[];
    combos_s=nchoosek(1:size(X{1},1),2);
    for i=1:length(X)
        v=[];
        net=X{i};
        [M,Q]=community_louvain(net);
        for ii=1:length(combos_s)
            m1=M(combos_s(ii,1),:);
            m2=M(combos_s(ii,2),:);
            if isequal(m1,m2)
                v=[v;1];
            else
                v=[v;0];
            end
        end
        Vs=cat(3,Vs,squareform(v));
    end
    Vs=sum(Vs,3);
    
    [M,Q]=community_louvain(sum(Vs,3));
    Opt_rank=max(M);
    
elseif Type==2
    Vs=[];
    combos_s=nchoosek(1:size(X{1},1),2);
    X=mat2tiles(X,[1,length(X)/Multi_type]);
        
    for it=1:combos_s
        c=[];
        net=X{it}';
        
        N=length(net{1});
        T=length(net);
        ii=[]; jj=[]; vv=[];
        twomu=0;
        clear AA all2all
        for s=1:T
            indx=[1:N]'+(s-1)*N;
            [i,j,v]=find(net{s});
            ii=[ii;indx(i)]; jj=[jj;indx(j)];vv=[vv;v];
            k=sum(net{s});
            kv=zeros(N*T,1);
            twom=sum(k);
            twomu=twomu+twom;
            kv(indx)=k/twom;
            kcell{s}=kv;
        end
        gamma=1;
        omega=1;
        AA = sparse(ii,jj,vv,N*T,N*T);
        clear ii jj vv
        kvec = full(sum(AA));
        all2all = N*[(-T+1):-1,1:(T-1)];
        AA = AA + omega*spdiags(ones(N*T,2*T-2),all2all,N*T,N*T);
        twomu=twomu+T*omega*N*(T-1);
        B = @(i) AA(:,i) - gamma*kcell{ceil(i/(N+eps))}*kvec(i);
        [S,Q] = iterated_genlouvain(B);
        
        S = reshape(S,N,T);
        S=S(:,1);
        c=[];
        for itv=1:length(combos_s)
            m1=S(combos_s(itv,1),:);
            m2=S(combos_s(itv,2),:);
            if isequal(m1,m2)
                c=[c;1];
            else
                c=[c;0];
            end
        end
        Vs=cat(3,Vs,squareform(c));
        
        
    end
       
    N=length(Vs(:,:,1));
    T=size(Vs,3);
    ii=[]; jj=[]; vv=[];
    twomu=0;
    for s=1:T
        indx=[1:N]'+(s-1)*N;
        [i,j,v]=find(Vs(:,:,s));
        ii=[ii;indx(i)]; jj=[jj;indx(j)]; vv=[vv;v];
        k=sum(Vs(:,:,s));
        kv=zeros(N*T,1);
        twom=sum(k);
        twomu=twomu+twom;
        kv(indx)=k/twom;
        kcell{s}=kv;
    end
    gamma=1;
    omega=1;
    AA = sparse(ii,jj,vv,N*T,N*T);
    clear ii jj vv
    kvec = full(sum(AA));
    all2all = N*[(-T+1):-1,1:(T-1)];
    AA = AA + omega*spdiags(ones(N*T,2*T-2),all2all,N*T,N*T);
    twomu=twomu+T*omega*N*(T-1);
    B = @(i) AA(:,i) - gamma*kcell{ceil(i/(N+eps))}*kvec(i);
    
    [S,Q] = iterated_genlouvain(B);
    S = reshape(S,N,T);
    Opt_rank=max(max(S));
end

end