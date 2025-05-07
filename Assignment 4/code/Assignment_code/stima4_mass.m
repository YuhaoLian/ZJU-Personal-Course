function M = stima4_mass(vertices)
vol = abs(det([ones(1,4); vertices'])) / 6;
M = (vol / 20) * (ones(4) + eye(4));
end