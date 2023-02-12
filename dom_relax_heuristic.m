function [obj,t,x, count] = dom_relax_heuristic(G, seed, varargin)
tic;
cplexpath;
rng(seed);

n = length(G);
f = ones(n,1);

A = [-1*(G+eye(n)); -1*eye(n)];
b = [-1*ones(n,1); zeros(n,1)];

Aeq = zeros(1,n);
beq = 0;

x=zeros(n,1);
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
while min(G*floor(x)+floor(x)) == 0
    
    if max(x) ~= 0
       
        k = x;
        k(fixed) = 0;
        k2 = k.^2;
        
        v = weightrand(k2);
        
        fixed = [fixed, v];
        
        Aeq(size(Aeq,1)+1,v) = 1;
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
obj = sum(x);
t = toc;
end