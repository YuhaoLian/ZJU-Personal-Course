% 精确解函数
function val = exact_solution(points)
    x = points(:,1); y = points(:,2);
    val = exp(2*x+2*y) .* sin(2*x+2*y);
end