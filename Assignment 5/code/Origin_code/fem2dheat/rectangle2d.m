function rectangle2d(n)
%Generate the mesh data file for a rectangle
%n=10;
xmin=-1;xmax=1; ymin=0;ymax=1;
nx=2*n-1; ny=n;
x=linspace(xmin,xmax,nx);
y=linspace(ymin,ymax,ny);
np=0; % no of vertices
X=zeros(2*nx*ny,1);
Y=zeros(2*nx*ny,1);
for k=1:ny
   for j=1:nx
       np  = np + 1;
       X(np) = x(j);
       Y(np) = y(k);
   end
end
for k=1:ny-1
   for j=1:nx-1
       np  = np + 1;
      X(np) = (x(j)+x(j+1))/2;
      Y(np) = (y(k)+y(k+1))/2;
   end
end
tri=delaunayTriangulation([X(1:np),Y(1:np)]);
%triplot(tri)
e=tri.edges;           % all edges
eb=tri.freeBoundary;   % boundary edges
nt   = size(tri,1); % no of triangles
fprintf(1,'Number of points         = %d\n', np);
fprintf(1,'Number of triangles      = %d\n', nt);

n_e  = size(e ,1); % no of all edges
n_eb = size(eb,1); % no of boundary edges
ebd = zeros(size(eb));
ebn = zeros(size(eb));
n_ebd = 0;
n_ebn = 0;
for j =1:n_eb
    yp1 = Y(eb(j,1));
    yp2 = Y(eb(j,2));
    if(yp1==1 && yp2==1)
        n_ebn = n_ebn+1;
        ebn(n_ebn,1) = eb(j,1);
        ebn(n_ebn,2) = eb(j,2);
    else
        n_ebd = n_ebd+1;
        ebd(n_ebd,1) = eb(j,1);
        ebd(n_ebd,2) = eb(j,2);
    end
end

ebn = ebn(1:n_ebn,:);
ebd = ebd(1:n_ebd,:);

fid=fopen('coordinates.dat','w');
fprintf(fid,'%24.14e %24.14e\n', tri.Points');
fclose(fid);

fid=fopen('elements3.dat','w');
fprintf(fid,'%8d %8d %8d\n', tri.ConnectivityList');
fclose(fid);

fid=fopen('dirichlet.dat','w');
fprintf(fid,'%8d %8d\n', ebd');
fclose(fid);

fid=fopen('neumann.dat','w');
fprintf(fid,'%8d %8d\n', ebn');
fclose(fid);
