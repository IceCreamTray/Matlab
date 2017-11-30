function [magn] = magnetization (Min, xdim, ydim)
for a=2:xdim-1
    for b=2:ydim-1
        Magn(a,b) = Min(a,b)*Min(a+1,b)+Min(a,b)*Min(a-1,b)+Min(a,b)*Min(a,b-1)+Min(a,b)*Min(a,b+1);
    end
end
magn = sum(sum(Magn)');
end