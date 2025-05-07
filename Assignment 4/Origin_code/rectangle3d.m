function rectangle3d(n)
%Generate the mesh data file for a cuboid
%n=10;
xmin=0;xmax=1;
ymin=0;ymax=1;zmin=0;zmax=1;
nx=n; ny=n; nz=n;
x=linspace(xmin,xmax,nx);
y=linspace(ymin,ymax,ny);
z=linspace(zmin,zmax,nz);
np=0; % no. of vertices
X=zeros(2*nx*ny*nz,1);
Y=zeros(2*nx*ny*nz,1);
Z=zeros(2*nx*ny*nz,1);
for l=1:nz
    for k=1:ny
        for j=1:nx
            np = np+ 1;
            X(np) = x(j);
            Y(np) = y(k);
            Z(np) = z(l);
        end
    end
end
for l=1:nz-1
    for k=1:ny-1
        for j=1:nx-1
            np = np+ 1;
            X(np) = (x(j)+x(j+1))/2;
            Y(np) = (y(k)+y(k+1))/2;
            Z(np) = (z(l)+z(l+1))/2;
        end
    end
end
tri=delaunayTriangulation([X(1:np),Y(1:np),Z(1:np)]);
coordinates = tri.Points; % all coordinate points
%element4 = tri.ConnectivityList; % all tetrahedrons
%e=tri.edges;
% all edges
eb=tri.freeBoundary;   % boundary faces
%tetramesh(tri,'FaceAlpha',0);

%%  绘图
% figure();
% trisurf(eb,coordinates(:,1),coordinates(:,2),coordinates(:,3),1);

%%
nt = size(tri,1); % no of triangles
fprintf(1,'Number of points = %d\n', np);
fprintf(1,'Number of triangles = %d\n', nt);

n_eb = size(eb,1); % no of boundary facets
ebd = zeros(size(eb));
ebn = zeros(size(eb));
n_ebd = 0;
n_ebn = 0;
tp = [1;1;1];
for j =1:n_eb
    t1 = sum(abs(X(eb(j,:))-tp));
    t2 = sum(abs(Y(eb(j,:))-tp));
    t3 = sum(abs(Z(eb(j,:))-tp));
    if(t1==0 || t2==0 ||t3 ==0)
        n_ebn = n_ebn+1;
        ebn(n_ebn,:) = eb(j,:);
    else
        n_ebd = n_ebd+1;
        ebd(n_ebd,:) = eb(j,:);
    end
end
ebn = ebn(1:n_ebn,:).';
ebd = ebd(1:n_ebd,:).';

fid=fopen('coordinates.dat','w');
fprintf(fid,'%24.14e %24.14e %24.14e\n',tri.Points');
fclose(fid);
fid=fopen('elements4.dat','w');
fprintf(fid,'%8d %8d %8d %8d\n',tri.ConnectivityList');
fclose(fid);
fid=fopen('dirichlet.dat','w');
fprintf(fid,'%8d %8d %8d\n', ebd);
fclose(fid);
fid=fopen('neumann.dat','w');
fprintf(fid,'%8d %8d %8d\n', ebn);
fclose(fid);