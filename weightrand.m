function index = weightrand(v)

sv = sum(v);
lv = length(v);
rn = rand*sv;
cv = cumsum(v);

for i=1:lv
    if(rn <= cv(i))
        index = i;
        return;
    end
end
