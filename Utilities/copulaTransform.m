function x = copulaTransform(x)
[r,c] = size(x);
disp('Copula Transform...')
for t=1:c,
[f,x1] = ecdf(x(:,t));
for t1=1:r,
   x(t1,t)=max(f(x1<=x(t1,t)));
end
end
end
