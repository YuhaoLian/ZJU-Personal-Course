function val = f(points)
    x = points(:,1); y = points(:,2);
    val = exp(2*x+2*y) .* (4*sin(2*x+2*y) - 16*cos(2*x+2*y));
end