# Polynomial Chaos Toolbox in Julia (Based on SymPy)
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

function innerProduct(P,Z,rangeZ)
    intP = [];
    for p in P
        for (i,z) in enumerate(Z)
            zmin = rangeZ[i,1];
            zmax = rangeZ[i,2];
            p = integrate(p,(z,zmin,zmax));
        end
        push!(intP,p);
    end
    return(reshape(intP,size(P)));
end

function getBasis(name::String,z::Array{SymPy.Sym,1},N::Integer)::Array{SymPy.Sym,1}
    nvar = length(z);
    @vars x
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

function symVec(z::String,nvar::Integer)
    Z = Array{SymPy.Sym,1}(undef,nvar);
    for i in 1:nvar
        Z[i] = symbols("$(z)$(i)",real=true);
    end
    return(Z);
end

Z = symVec("z",1);
PDF = 0.5^length(Z);
rangeZ = repeat([-1 1], length(Z));
Phi = getBasis("Legendre",Z,6);
intP1 = innerProduct(Phi*PDF,Z,rangeZ);

f = sin(Z[1]^3);
@time X = innerProduct(f*Phi*0.5,Z,rangeZ); # Really slow.
fhat = X'*Phi;
Xval = -1:0.01:1;
f1val = [subs(fhat,Z[1]=>x) for x in Xval];
fval =  [subs(f,Z[1]=>x) for x in Xval];

using Plots
plot(Xval,fval);
p = plot!(Xval,f1val);
plot(p,label=["f","fhat"])

@time intP1 = innerProduct(Phi.*Phi*PDF,Z,rangeZ);
@time W = innerProduct(Phi*Phi'*PDF,Z,rangeZ);
print("Done")
