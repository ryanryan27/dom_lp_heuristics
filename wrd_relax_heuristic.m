function [obj,t,xw, xww] = wrd_relax_heuristic(G,seed)


tic;
cplexpath;

rng(seed);


[A,b,Aeq,beq,f,lb,ub] = wrd_rlp_build(G);
n=length(G);
x=zeros(2*n,1);
known = zeros(2*n,1);



while sum(wrd_penalty(floor(x(1:n))+2*floor(x(n+1:2*n)),G,1,1)) > 0
    if max(x)~=0
        
       k = x;
       k(find(known)) = 0;
       k2 = k.^2;
       v = weightrand(k2);

        if v <= n
            v2 = v + n;
        else
            v2 = v - n;
        end

        known(v2) = 1;
        known(v)=1;
        Aeq(size(Aeq,1)+1,v) = 1;
        Aeq(size(Aeq,1),v2) = 1;
        beq(size(Aeq,1),1) = 1;

    end

    options = cplexoptimset;
    options.Algorithm = 'primal';
    options.Display = 'off';

    xy = cplexlp(f,A,b,Aeq,beq,lb,ub,[],options);

    x=xy(1:2*n);

end
xw=floor(x);
xww=xw(1:n)+ 2*xw(n+1:2*n);
obj2 = sum(x(1:n)) + 2*sum(x(n+1:2*n));
obj = sum(xww);
if(abs(obj2 - obj) > 1e-4)
    disp('STRANGE - OBJ AND OBJ2 ARE DIFFERENT?');
    pause
end
t = toc;

end