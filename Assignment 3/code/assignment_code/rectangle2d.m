function rectangle2d(n)
xmin = -1; xmax = 1; ymin = 0; ymax = 1;
nx = 2*n-1; ny = n;
x = linspace(xmin, xmax, nx);
y = linspace(ymin, ymax, ny);

% Generate nodes
np = 0;
X = zeros(2*nx*ny, 1); Y = zeros(2*nx*ny, 1);
for k = 1:ny
    for j = 1:nx
        np = np + 1;
        X(np) = x(j); Y(np) = y(k);
    end
end
for k = 1:ny-1
    for j = 1:nx-1
        np = np + 1;
        X(np) = (x(j)+x(j+1))/2;
        Y(np) = (y(k)+y(k+1))/2;
    end
end

% Triangulate and mark all boundaries as Dirichlet
tri = delaunayTriangulation(X(1:np), Y(1:np));
eb = tri.freeBoundary;

% Save mesh data
dlmwrite('coordinates.dat', tri.Points, 'precision', '%24.14e');
dlmwrite('elements3.dat', tri.ConnectivityList, 'delimiter', '\t');
dlmwrite('dirichlet.dat', eb, 'delimiter', '\t');
dlmwrite('neumann.dat', [], 'delimiter', '\t'); % Empty Neumann
end