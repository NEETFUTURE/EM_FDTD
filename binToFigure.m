function data = binToFigure(filename,Nx,Ny)
fid =  fopen(filename);
data = fread(fid, [Nx Ny], 'float');
figure
surf(data)
zlim([-0.1 0.1])
end