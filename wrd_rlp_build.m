function [A,b,Aeq,beq,f,lb,ub,ctype] = wrd_rlp_build(G)

n = length(G);
f = [ones(n,1); 2.*ones(n,1); zeros(n*n,1)];
lb = zeros(size(f));
ub = ones(size(f));
ctype = '';
for i=1:2*n
    ctype = [ctype 'C'];
end
for i=1:n^2
    ctype = [ctype 'C'];
end

m = sum(sum(G))/2;

count = 0;


twoapartcount = 0;
for i=1:n
    neighbours = find(G(i,:));
    twoindex = zeros(n,1);
    for j=1:length(neighbours)
        neighbours2 = find(G(neighbours(j),:));
        for k=1:length(neighbours2)
            if(neighbours2(k) ~= i && sum(neighbours2(k) == neighbours) == 0)
                twoindex(neighbours2(k)) = 1;
            end
        end
    end
    twoapartcount = twoapartcount + sum(twoindex);
end


A = sparse(n+2*m+twoapartcount,n*(n+2));
b = sparse(n+2*m+twoapartcount,1);

%constraint 1
A(1:n,:) = [-1*(G+eye(n)), -1*(G+eye(n)), zeros(n,n*n)];
b(1:n) = -1*ones(n,1);
count = n;

%constraint 2
for i = 1:n
    c = find(G(i,:));
    
    for j = c
       count = count + 1;
       A(count,i) = -1; %x vars
       A(count,i+n) = -1; %w vars
       A(count,(i+1)*n + j) = 1;
      
       b(count) = 0;
    end
end

%constraint 3
for i=1:n
    for j=two_apart(G,i)
        %D = zeros(1,n*(n+1));
        
        ni = find(G(i,:));
        
        cni = [i, ni];
        count = count + 1;
        A(count,cni) = -1; %x vars
        A(count,cni+n) = -2; %w vars
        for k=ni
            if(G(k,j))
                A(count,(k+1)*n+j) = 1;
            end
        end
        
        %A(count,:) = D;
        b(count) = -1;
    end
end

Aeq = sparse(n,n*(n+2));
beq = sparse(n,1);
counteq = 0;



%constraint 4
for j = 1:n
    d = find(G(j,:));
    %D = zeros(1,n*(n+1));
    counteq = counteq + 1;
    for i = d
        Aeq(counteq,(i+1)*n + j) = 1;
    end
    Aeq(counteq,j) = 1; %x vars
    Aeq(counteq,j+n) = 1; %w vars
    beq(counteq) = 1;

end

%A = A(1:count,:);
%b = b(1:count);

if count ~= n + 2*m + twoapartcount
    count
    n + 2*m + twoapartcount
    disp('wrong')
    pause;
end



end

