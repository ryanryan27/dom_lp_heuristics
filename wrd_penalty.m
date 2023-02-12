

function UG = wrd_penalty(x, G, sparse_d, sparse_i)

sparse = 1;
type = 0;

n = size(G,1);

if type
    
    UG = zeros(1,n);
    for i=1:n
        if(x(i) > 0)
            continue;
        end
        
        guarded = 0;
        if ~sparse
            for j=1:n
                if(G(i,j))
                    if(x(j) == 2)
                        guarded = 1;
                        break;
                    end
                    if(x(j) == 1)
                        
                        kguarded = 1;
                        for k=1:n
                            if(G(j,k) && x(k) == 0 && k ~= i && G(i,k)==0)
                                kguarded = 0;
                                for l=1:n
                                    if(G(k,l) && x(l) && l ~= j)
                                        kguarded = 1;
                                        break;
                                    end
                                end
                                if(kguarded == 0)
                                    break;
                                end
                            end
                        end
                        if(kguarded == 0)
                            continue;
                        else
                            guarded = 1;
                            break;
                        end
                    end
                end
            end
            if(guarded == 0)
                UG(i) = 1;
            end
        else
            ni = sparse_d(sparse_i(i):sparse_i(i+1)-1);
            for jj=1:length(ni)
                j = ni(jj);
                if(x(j) == 2)
                    guarded = 1;
                    break;
                end
                if(x(j) == 1)
                    
                    kguarded = 1;
                    nj = sparse_d(sparse_i(j):sparse_i(j+1)-1);
                    for kk=1:length(nj)
                        k=nj(kk);
                        if(x(k) == 0 && k ~= i && G(i,k)==0)
                            kguarded = 0;
                            nk = sparse_d(sparse_i(k):sparse_i(k+1)-1);
                            for ll = 1:length(nk)
                                l = nk(ll);
                                if(x(l) && l ~= j)
                                    kguarded = 1;
                                    break;
                                end
                            end
                            if(kguarded == 0)
                                break;
                            end
                        end
                    end
                    if(kguarded == 0)
                        continue;
                    else
                        guarded = 1;
                        break;
                    end
                end
                
            end
            if(guarded == 0)
                UG(i) = 1;
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    
    % vertices with exactly one guard
    g = find(x == 1);
    k = length(g);
    
    % initialise every vertex to be unguarded
    %%UG = ones(k + 2, length(x));
    UG = ones(1,length(x));
    
    % for each vertex with exactly one guard
    for ii = 1:k
        i = g(ii);
        
        % get  its unguarded neighbours
        
        ni = find((G(:,i) - x) == 1);
        
        
        
        % temporarily assume that i has no guard, i.e. it moves to some j
        x(i) = 0;
        
        % get the adjacencies of its unguarded neighbours
        GG = G(ni,:);
        
        % for each neighbour with no guard
        for jj = 1:length(ni)
            j = ni(jj);
            
            % the guard at i moved to j
            x(j) = 1;
            
            % since j has a guard now, this is used in the next bit to
            % ensure that it doesnt think j is unguarded
            GG(jj,j) = 1;
            
            % GG*x gives a vertex with each value corresponding to whether
            % or not a particular neighbour of i is still guarded after the
            % guard moves to j
            if min(GG*x) == 0
                % some neighbour of i is undefended after the guard moves
                %%UG(ii,j) = 1;
                
            else
                %%UG(ii,j) = 0;
                UG(j) = 0;
            end
            
            %undo the temp changes to GG and x
            GG(jj,j) = 0;            
            x(j) = 0;
        end
        
        %put the guard back at i
        x(i) = 1;
        
        
        
    end
    
    % find vertices with exactly two guards
    twos = find(x == 2);
    
    % for each vertex with two guards
    for ii = 1:length(twos)
        i = twos(ii);
        
        % every neighbour of i is guarded
        %UG(k+2,:) = UG(k+2,:).* ~G(i,:);
        UG = UG.*~G(i,:);
    end
    
    % every vertex with a guard is defended
    %%UG(k+1,:) = (x == 0)';
    %%UG = min(UG);
    UG = UG.* (x == 0)';
end
end
