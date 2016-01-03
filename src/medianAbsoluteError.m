% Downloaded from http://www.run.montefiore.ulg.ac.be/~liao/DMFSGD.html

function m = medianAbsoluteError(X,Xhat)

ids = find(X>0);

x = X(ids);
xhat = Xhat(ids);

m = median(abs(x-xhat));

return