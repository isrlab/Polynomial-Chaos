# Polynomial Chaos Toolbox in Julia 1.0
# (c) Raktim Bhattacharya
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# Contributors:
# 1. Raktim Bhattacharya -- raktim@tamu.edu -- Aerospace Engineering, Texas A&M University
#
#
#
#
################################################################################

using DynamicPolynomials, Combinatorics

function legendrePolynomials(x::PolyVar,N::Integer)::Array{Polynomial,1}
    P =  Array{DynamicPolynomials.Polynomial}(N+1);
    for n in 0:N
        if n == 0
            p = 0.0+x^0;
        elseif n == 1
            p = 0.0+x;
        else
            m = n-1;
            p = ((2*m+1)*x*P[n] -m*P[n-1])/n;
        end
        P[n+1] = p;
    end
    return P
end

function integratePolynomial(p::Polynomial,x::PolyVar)::Polynomial
    newP = Polynomial{true,Float64}(p);
    i = find(isequal(x),newP.x.vars);
    for j = 1:length(newP.x.Z)
        newPow = newP.x.Z[j][i][1]+1;
        newP.x.Z[j][i] = newPow;
        newP.a[j] = newP.a[j]/newPow;
    end
    return(newP);
end

function integratePolynomial(m::Monomial,x::PolyVar)::Polynomial
    return(integratePolynomial(Polynomial(m),x))
end

function integratePolynomial(t::Term,x::PolyVar)::Polynomial
    return(integratePolynomial(Polynomial(t),x))
end

function polyVal(p::Polynomial,z::PolyVar,zVal::Float64)::Polynomial
# Need to implement fast polynomial evaluation routine.
# Subs is too slow.
end

# Integrate polyvariate polynomials in z over zRange.
function innerProduct(P,Z,rangeZ)
    Q = [];

    for p in P
        for z in Z
            p = integratePolynomial(p,z);
        end
        push!(Q,p)
    end

    # Now compute definite integral, sub numerical values from rangeZ.
    for (k,q) in enumerate(Q)
        for (i,z) in enumerate(Z)
            Q[k] = subs(Q[k],z=>rangeZ[i,2]) - subs(Q[k],z=>rangeZ[i,1]); # This could be slow.
        end
    end

    Q = convert(Array{Real,1},Q);
    Q = reshape(Q,size(P));
    return(Q);
end

function getBasis(name::String,z::Array{DynamicPolynomials.PolyVar{true},1},N::Integer)::Array{DynamicPolynomials.Polynomial{true,Float64},1}
    nvar = length(z);
    @polyvar x
    if name=="Legendre"
        basis = legendrePolynomials(x,N); # Add other basis functions
    end

    Phi = [];
    for n in 0:N
        for indices in multiexponents(nvar,n)
            polyBasis = 1;
            for (i,index) in enumerate(indices)
                polyBasis = polyBasis*subs(basis[index+1],x=>z[i]);
            end
            push!(Phi,polyBasis);
        end
    end
    return(Phi);
end

@polyvar z[1:2];

Phi = getBasis("Legendre",z,10);
rangeZ = repmat([-1 1],length(z));

PDF = (0.5)^length(z);
intQ1 = innerProduct(Phi*PDF,z,rangeZ)
