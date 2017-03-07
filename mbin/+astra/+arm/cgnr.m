function x = cgnr(A,b,niter,x0)

x = x0;
r = b - A*x;
z = A'*r;
p = z;

for i = 1:niter
    w     = A*p;
    eta   = z'*z;
    alpha = eta / (w'*w);
    x     = x + alpha*p;
    r     = r - alpha*w;
    z     = A'*r;
    kappa = z'*z;
    beta  = kappa/eta;
    p     = z + beta*p;
end