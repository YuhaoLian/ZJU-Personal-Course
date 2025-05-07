function M = stima4(vertices)
% Given the vertices of a tetrahedron, evaluate the stiffness matrix element
G = [ones(1,4);vertices'] \ [zeros(1,3);eye(3)];
M = det([ones(1,4);vertices']) * G * G'/6;
end
