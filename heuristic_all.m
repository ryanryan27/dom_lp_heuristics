function [obj,x,t] = heuristic_all(G,dom_type, iterations)
tstart = tic;
n = length(G);
obj = 2*n+1;
x = zeros(n,1);

for i = 1:iterations

    
    switch dom_type   
        case 'd'
            %domination
            [obji,~,xi] = dom_relax_heuristic(G,i);
        case 't'
            %total domination
            [obji,~,xi] = td_relax_heuristic(G,i);
        case 's'
            %secure domination
            [obji,~,xi] = sd_relax_heuristic(G,i);
        case 'r'
            %roman domination
            [obji,~,~,xi] = rd_relax_heuristic(G,i);
        case 'w'
            %weak roman domination
            [obji,~,~,xi] = wrd_relax_heuristic(G, i);
        otherwise
            %default option is domination
            [obji,~,xi] = dom_relax_heuristic(G,i);
    end
    if obji < obj
        obj = obji;
        x = xi;
    end
end
t = toc(tstart);
end