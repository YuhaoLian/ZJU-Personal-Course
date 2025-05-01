function M = stima3 ( vertices )
% Given the vertices of a triangle, evaluate the
% stiffness matrix element
G = [ ones(1,3); vertices' ] \ [ zeros(1,2); eye(2) ];
M = det ( [ ones(1,3); vertices' ] ) * G * G' / 2;
end