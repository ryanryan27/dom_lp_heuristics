function [obj,t,x,xw, count] = rd_relax_heuristic(G, seed, varargin)
tic;
cplexpath;
rng(seed);

n = length(G);
f = [ones(n,1); 2*ones(n,1)];

A = [-1*(eye(n)),-1*(G + eye(n)); -1*eye(2*n)];
b = [-1*ones(n,1); zeros(2*n,1)];

Aeq = zeros(1,2*n);
beq = 0;

x=zeros(2*n,1);
fixed = [];

count = 0;

if length(varargin) > 0
   fixed = varargin{1}; 
end

%x(fixed,1) = 1;

for i = 1:length(fixed)
    j = size(Aeq,1)+1;
    Aeq(j,i) = 1;
    beq(j,1) = 1;
end



while min((G)*floor(x((n+1):2*n))+floor(x(1:n)) + floor(x((n+1):2*n))) == 0
    
    if max(x) ~= 0
       
        k = x;
        k(fixed) = 0;
        k2 = k.^2;

        v = weightrand(k2);

        if v <= n
            v2 = v + n;
        else
            v2 = v - n;
        end
            
        fixed = [fixed, v, v2];
        Aeq(size(Aeq,1)+1,v) = 1;
        Aeq(size(Aeq,1),v2) = 1;
        beq(size(Aeq,1),1) = 1;

            
    end
    
    options = cplexoptimset;
    options.Algorithm = 'primal';
    options.Display = 'off';
        
    x = cplexlp(f,A,b,Aeq,beq,[],[],[],options);
    
    xx = floor(x);

    
    count = count + 1;
    
end
x = floor(x);
xw = x(1:n) + 2*x((n+1):2*n);
obj = sum(x(1:n)) + 2*sum(x((n+1):2*n));
t = toc;
end