% Downloaded from http://www.run.montefiore.ulg.ac.be/~liao/DMFSGD.html

function  s = stress(X,Xa,isL2)

if ~exist('isL2')
    isL2 = 1;
end

ids = find(X>1e-5);

x = X(ids);
xa = Xa(ids);

d = x-xa;
if isL2>0
    s = sqrt((d'*d)/((x'*x)+1e-8));
else
    s = sum(abs(d))/(sum(x)+1e-8);
end
return