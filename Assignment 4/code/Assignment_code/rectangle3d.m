function rectangle3d(n)
xmin = -1; xmax = 1;
ymin = -1; ymax = 1;
zmin = -1; zmax = 1;
nx = n; ny = n; nz = n;
x = linspace(xmin, xmax, nx);
y = linspace(ymin, ymax, ny);
z = linspace(zmin, zmax, nz);
np = 0;
X = zeros(2*nx*ny*nz, 1);
Y = zeros(2*nx*ny*nz, 1);
Z = zeros(2*nx*ny*nz, 1);

% Generate vertex points
for l = 1:nz
    for k = 1:ny
        for j = 1:nx
            np = np + 1;
            X(np) = x(j);
            Y(np) = y(k);
            Z(np) = z(l);
        end
    end
end

% Generate midpoints
for l = 1:nz-1
    for k = 1:ny-1
        for j = 1:nx-1
            np = np + 1;
            X(np) = (x(j) + x(j+1))/2;
            Y(np) = (y(k) + y(k+1))/2;
            Z(np) = (z(l) + z(l+1))/2;
        end
    end
end

tri = delaunayTriangulation([X(1:np), Y(1:np), Z(1:np)]);
coordinates = tri.Points;
eb = tri.freeBoundary;

% All boundaries are Dirichlet
dirichlet = eb;
neumann = [];

% Save mesh files
fid = fopen('coordinates.dat', 'w');
fprintf(fid, '%24.14e %24.14e %24.14e\n', coordinates');
fclose(fid);

fid = fopen('elements4.dat', 'w');
fprintf(fid, '%8d %8d %8d %8d\n', tri.ConnectivityList');
fclose(fid);

fid = fopen('dirichlet.dat', 'w');
fprintf(fid, '%8d %8d %8d\n', dirichlet');
fclose(fid);

fid = fopen('neumann.dat', 'w');
fprintf(fid, '%8d %8d %8d\n', neumann');
fclose(fid);
end