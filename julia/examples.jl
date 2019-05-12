# Code example to use the polynomial chaos toolbox.
# Function approximation

include("./PolynomialChaos.jl");

# Example 1
Z = symVec("z",2);
PDF = 0.5^length(Z);
rangeZ = repeat([-1 1], length(Z));
Phi = getBasis("Legendre",Z,5);

@time M1_Phi = innerProduct(Phi*PDF,Z,rangeZ);
@time M2_Phi = innerProduct(Phi.*Phi*PDF,Z,rangeZ);
@time W = innerProduct(Phi*Phi'*PDF,Z,rangeZ);
println("Done")


# Example 2 -- More complicated functions
Z = symVec("z",1);
PDF = 0.5^length(Z);
rangeZ = repeat([-1 1], length(Z));
Phi = getBasis("Legendre",Z,3);
f = sin(Z[1]^3);
@time X = innerProduct(f*Phi*0.5,Z,rangeZ); # Really slow -- try non intrusive.
fhat = X'*Phi;
Xval = -1:0.01:1;
f1val = [subs(fhat,Z[1]=>x) for x in Xval];
fval =  [subs(f,Z[1]=>x) for x in Xval];

using PyPlot
plot(Xval,fval);
p = plot(Xval,f1val);
legend(["f","fhat"]);
