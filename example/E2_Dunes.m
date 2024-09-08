% Load the training image
load('E2 - Dunes/TI_dunes.mat')

%% Image size
[nx, ny] = size(ti);

%Display the image
figure(8)
imagesc(ti)
colorbar()
colormap 'gray'
clim([1 3])
set(gca,'YDir','normal')

% Reshape to a vector
Z = reshape(ti,[],1);
%% Generate a 2D grid
x0 = grille2(1,nx,1,1,ny,1);

%% Compute the transiogram, for categorical data
categ = 1;
display = 1;
icode = 6;
[gh, nh]=GeoStatFFT(x0,Z,icode,categ,display);

figure(2)
for i=1:3
    for j=1:3
        subplot(3,3,((i-1)*3+j))
        plot(gh{i,j}(nx:end,ny),'-k')
    end
end
%% Bivariate probabilities, for categorical data
categ = 1;
display = 1;
icode = 5;
[gh, nh]=GeoStatFFT(x0,Z,icode,categ,display);

figure(2)
for i=1:3
    for j=1:3
        subplot(3,3,((i-1)*3+j))
        plot(gh{i,j}(nx:end,ny),'-k')
    end
end

%% Non-ergodic transiograms
categ = 1;
display = 1;
icode = 7;
[gh, nh]=GeoStatFFT(x0,Z,icode,categ,display);

figure(2)
for i=1:3
    for j=1:3
        subplot(3,3,((i-1)*3+j))
        plot(gh{i,j}(nx:end,ny),'-k')
    end
end

%% TOC 
categ = 1;
display = 1;
icode = [11 5 0];
figure(1)
[gh, nh]=GeoStatFFT(x0,Z,icode,categ,display);




icode = [11 5 5];
figure(2)
[gh, nh]=GeoStatFFT(x0,Z,icode,categ,display);



icode = [11 0 5];
figure(3)
[gh, nh]=GeoStatFFT(x0,Z,icode,categ,display);



icode = [11 -5 5];
figure(4)
[gh, nh]=GeoStatFFT(x0,Z,icode,categ,display);



