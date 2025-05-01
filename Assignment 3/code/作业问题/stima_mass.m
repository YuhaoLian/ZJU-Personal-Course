function M = stima_mass(vertices)
area = det([1 1 1; vertices']) / 2; % Triangle area
M = (area / 12) * [2 1 1; 1 2 1; 1 1 2]; % Mass matrix
end