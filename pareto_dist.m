function result = pareto_dist(xmin, alpha)
    result = ceil(xmin/(1-rand)^(1/alpha)-xmin);
