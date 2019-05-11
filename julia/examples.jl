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

using PyPlot
plot(Xval,fval);
p = plot(Xval,f1val);
legend(["f","fhat"]);

@time intP1 = innerProduct(Phi.*Phi*PDF,Z,rangeZ);
@time W = innerProduct(Phi*Phi'*PDF,Z,rangeZ);
print("Done")
