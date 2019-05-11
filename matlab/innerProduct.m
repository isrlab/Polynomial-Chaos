function val = innerProduct(fcn,vars,lim)

val = fcn;
for i=1:length(vars)
    val = int(val,vars(i),lim(i,1),lim(i,2));
end

    

