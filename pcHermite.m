function P = pcHermite(x,n)

f1 = exp(-x^2/2);
f2 = exp(x^2/2);
f3 = (-1)^n;

P = expand(f3*simplify(f2*diff(f1,n)));