function M = stima3(vertices)
G = [ones(1,3); vertices'] \ [zeros(1,2); eye(2)];
M = det([ones(1,3); vertices']) * G * G' / 2; % Laplacian stiffness
end