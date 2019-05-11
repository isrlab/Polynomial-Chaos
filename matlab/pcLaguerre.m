function P = pcLaguerre(x,n)
if n>21 
    disp('Accuracy guaranteed for n<=21 (factorial)');
end

f = exp(-x)*x^n;
P = simplify(exp(x)*diff(f,n))/factorial(n);
P = expand(P);

