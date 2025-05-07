function M = stima4_mass(vertices)
d = det([ones(1,4); vertices']);
V = abs(d) / 6;
M = V * [2 1 1 1; 1 2 1 1; 1 1 2 1; 1 1 1 2] / 20;
end