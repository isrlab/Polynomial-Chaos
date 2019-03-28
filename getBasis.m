function [phi,phiComponents] = getBasis(type,n,p,varname);
% SYNTAX
% function [phi] = getBasis(type,n,p,varname);
% type = 'Hermite', 'Legendre', 'Jacobi'
% n = number of random variables
% p = maximum order of polynomials
% will give factorial(n+p)/(factorial(n)*factorial(p)) vector, phi
% and weighting function, W

for kk=1:n,
    x(kk,1)=sym([varname num2str(kk)],'real');
end
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

switch type
    case 'Hermite'
        W=1/sqrt((2*pi)^n)*exp(-(x'*x)/2);
        lowlim=-inf;uplim=inf;
    case 'Jacobi'
        disp(['This function is not implemented yet ' num2str(i)]);
    case 'Legendre'
        W=(1/2)^n;
        lowlim=-1;uplim=1;
    otherwise
        %type='Legendre';
        %W=(1/2)^n;
        error(['Unknown basis function:' type]);
end

% Now find all orthogonal polynomials from order 0 to p.
Pnum=0;
for i=0:p
%    disp(' ');
%     switch type
%         case 'Hermite'
%             disp(['Multivariate Hermite Polynomials of order ' num2str(i)]);
%         case 'Jacobi'
%             disp(['This function is not implemented yet ' num2str(i)]);
%         case 'Legendre'                
%             disp(['Multivariate Legendre Polynomials of order ' num2str(i)]);
%     end        
    ii = find(S==i);
    for j=1:length(ii)
        Poly = sym('1');
        for k=1:n
            switch type
                case 'Hermite'
                    Poly = Poly*pcHermite(x(k),M(ii(j),k));  
                case 'Jacobi'
                    Poly = Poly*pcJacobi(x(k),M(ii(j),k));
                case 'Legendre'
                    Poly = Poly*pcLegendre(x(k),M(ii(j),k));
            end
        end
        %pretty(simplify(expand(Poly)));
        Pnum=Pnum+1;
        %keyboard;
        phi(Pnum,1)=Poly;
    end
end
