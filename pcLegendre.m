function P = pcLegendre(x,n)
if n>21
    disp('Accuracy guaranteed for n<=21 (factorial)');
end

if n==0
    P = sym('1');
else
    v1 = 2^n;
    v2 = factorial(n);
    f = (x^2-1)^n;
    P = (diff(f,n)/v1)/v2;
    P = expand(P);
end

