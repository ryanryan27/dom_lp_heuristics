function [obj,x, count] = sd_relax_heuristic(G,seed, varargin)
 version = 0;
if length(varargin) > 0
    version = varargin{1};
end

cplexpath;
rng(seed);
epsilon = 0.0001;

[A,b,Aeq,beq,f,lb,ub] = sd_rlp_build(G);

n=length(G);
x=zeros(n,1);
known = zeros(n,1);
type = 1;
count = 0;
while sum(wrd_penalty(floor(x),G,1,1)) > 0
    if max(x)~=0
        if(type == 1)
           k = x;
           k(find(known)) = 0;
           k2 = k.^2;
           v = weightrand(k2);
        elseif(type == 2)
            k = x - 2*known;
            
            if version == 1
                z = find(k==max(k));
            else
                z = find(abs(k-max(k)) < epsilon);
            end
            
            count = count + 1;
            %fprintf([num2str(count) ' ' num2str(max(k)) ' ' num2str(length(z)) '\n']);
           v = z(randi(length(z)));
           
        end
      
        
        %v=z(1);
        known(v)=1;
        Aeq(size(Aeq,1)+1,v) = 1;
        beq(size(Aeq,1),1) = 1;
    end
        
    %much faster with milp for some reason?
    options = cplexoptimset;
    options.Algorithm = 'primal';
    options.Display = 'off';
    %xy = cplexmilp(f,A,b,Aeq,beq,[],[],[],lb,ub,ctype);
    
    xy = cplexlp(f,A,b,Aeq,beq,lb,ub,[],options);
    
    x=xy(1:n);
    
    count = count + 1;
end
x=floor(x);
obj = sum(x);
end