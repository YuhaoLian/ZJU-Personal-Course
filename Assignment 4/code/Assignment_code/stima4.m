function M = stima4(vertices)
d = det([ones(1,4); vertices']);
G = [ones(1,4); vertices'] \ [zeros(1,3); eye(3)];
M = abs(d) * G * G' / 6;
end