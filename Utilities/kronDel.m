function d = kronDel(y)
[R,~] = size(y);
d = zeros(R,R);
for r=1:R,
    for c=1:R,
        if y(r) == y(c)
            d(r,c) = 1;
        else
            d(r,c) = 0;
        end
    end
end
