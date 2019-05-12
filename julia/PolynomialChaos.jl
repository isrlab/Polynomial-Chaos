using SymPy, Combinatorics, LinearAlgebra

function legendrePolynomials(x::SymPy.Sym,N::Integer)::Array{SymPy.Sym,1}
    P =  Array{SymPy.Sym}(undef,N+1);
    for n in 0:N
        if n == 0
            p = 1;
        elseif n == 1
            p = x;
        else
            m = n-1;
            p = ((2*m+1)*x*P[n] -m*P[n-1])/n;
        end
        P[n+1] = p;
    end
    return P
end


# Have to make it typed.
function innerProduct(P,Z,rangeZ)
    intP = zeros(prod(size(P)),1);
    for k=1:length(P)
        p = P[k];
        for i=1:length(Z);
            z = Z[i];
            zmin = rangeZ[i,1];
            zmax = rangeZ[i,2];
            p = integrate(p,(z,zmin,zmax));
        end
        intP[k] = p; #Float64(p);
    end
    X = reshape(intP,size(P));
    return(X);
end

function getBasis(name::String,z::Array{SymPy.Sym,1},N::Integer)::Array{SymPy.Sym,1}
    nvar = length(z);
    @vars x
    if name=="Legendre"
        basis = legendrePolynomials(x,N); # Add other basis functions
    else
        error("Not implemented");
    end

    Phi = [];
    for n in 0:N
        for indices in multiexponents(nvar,n)
            polyBasis = 1;
            for (i,index) in enumerate(indices)
                polyBasis = simplify(polyBasis*subs(basis[index+1],x=>z[i]));
            end
            push!(Phi,polyBasis);
        end
    end
    return(Phi);
end

function symVec(z::String,nvar::Integer)::Array{SymPy.Sym,1}
    Z = Array{SymPy.Sym,1}(undef,nvar);
    for i in 1:nvar
        Z[i] = symbols("$(z)$(i)",real=true);
    end
    return(Z);
end
