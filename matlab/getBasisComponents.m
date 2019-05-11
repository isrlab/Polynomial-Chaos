function phiComponents = getBasisComponents(n,p)
% SYNTAX
% function [phi] = getBasis(type,n,p,varname);
% type = 'Hermite', 'Legendre', 'Jacobi'
% n = number of random variables
% p = maximum order of polynomials
% will give factorial(n+p)/(factorial(n)*factorial(p)) vector, phi
% and weighting function, W

endIndex = (p+1)^n-1;
v1 = [0:endIndex];
for i=1:length(v1)
    v2 = v1(i);
    for j=(n-1):(-1):0
        M(i,j+1) = floor(v2/(p+1)^j);
        v2 = mod(v2,(p+1)^j);
    end
end
S = sum(M,2);

Pnum=0;
for i=0:p
    ii = find(S==i);
    for j=1:length(ii)
        componentPoly = ones(1,n);
        for k=1:n
            componentPoly(1,k) = M(ii(j),k)+1;
        end
        Pnum=Pnum+1;
        phiComponents(Pnum,:) = componentPoly;
    end    
end

