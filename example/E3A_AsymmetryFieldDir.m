
%% Clear temporary variables
clear variables
% Load and process the png image. We don't have the real value so we
% extract it from the RBG code.

X = imread('E3 - SyntheticAsymmetry/Directional_field_v1.png');
[nx,ny, dim] = size(X );
x0 = grille2(1,nx,1,1,ny,1);
X = reshape(im2double(X(:,:,1:3)), [],3);

figure(1)
imagesc(reshape(X,[nx ny 3]));
colorbar();
colormap 'jet'
clim([-3.2 3.2])
set(gca,'YDir','normal')

[~,idx] = sortrows(rgb2hsv(X), -1);
C = gray(nx*ny);
C(idx,:) = C;

Im = C(:,1);
Im(Im<=0.04) = nan;
Im(Im>=0.96) = nan;
Im = reshape(Im,[nx ny]);

filledImage = fillmissing(Im, 'nearest'); % Fill NaNs with nearest neighbor

Z =  reshape(filledImage,[nx*ny 1]);
Z = (Z - mean(Z))/std(Z) + rand(nx*ny,1)/10;
A = reshape(Z,[nx ny]);
B = [1 1 1; 1 1 1; 1 1 1]/9;
Z = ECDF(reshape(conv2(A, B, "same"),[nx*ny 1])) ;

figure(2)
imagesc(reshape(Z(:,1),[nx ny]));
colorbar();
colormap 'jet'
clim([0 1])
set(gca,'YDir','normal')
fontsize(gca, 12, 'points')
%% Compute the directional asymetry
categ = 0;
display = 1;
rank = 1;
icode = 8;
[gh, nh] = GeoStatFFT(x0, Z, icode, categ, display, rank);

% Simulation of random fields
[z]=fftma([6 50 50 50 0 0 0],1,45124241,100,nx,1,ny,1);
for i =1:100
    gh1{i} = GeoStatFFT(x0, z(:,i), icode, categ, 0, rank);
end

figure(3)
for i=1:100
    plot(0:length(gh1{i}{1,1}(nx:end,ny))-1, gh1{i}{1,1}(nx:end,ny),'color', [0.9 0.9 0.9], LineWidth=1)
    hold on
    gmean(:,i) = gh1{i}{1,1}(nx:end,ny);
end
plot(0:length(mean(gmean,2))-1,mean(gmean,2), '-k',LineWidth=1.5)
hold on
plot(0:length(quantile(gmean,0.05,2))-1,quantile(gmean,0.05,2),'--k',LineWidth=1)
hold on
plot(0:length(quantile(gmean,0.95,2))-1,quantile(gmean,0.95,2),'--k',LineWidth=1)
hold on
h1 = plot(0:length(gh{1,1}(nx:end,ny))-1,gh{1,1}(nx:end,ny),'-g',LineWidth=2);
hold on
h2 = plot(0:length(gh{1,1}(nx,ny:end))-1,gh{1,1}(nx,ny:end),'-r',LineWidth=2);
legend([h1 h2],{'Vertical', 'Horizontal'})
xlabel('Distance')
ylabel('Directional asymmetry value')
fontsize(gca, 14, 'points')
xlim([0 floor(nx/4)])
ylim([-0.0495 0.0495])
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')


%% Compute the bivariate asymetry
categ = 0;
display = 1;
icode = 9;
[gh, nh]=GeoStatFFT(x0,Z,icode,categ,display,1);

% Simulation of random fields
[z]=fftma([6 50 50 50 0 0 0],1,45124241,100,nx,1,ny,1);
for i =1:100
    gh1{i} = GeoStatFFT(x0,z(:,i),icode,categ,0,1);
end

figure(4)
for i=1:100
    plot(0:length(gh1{i}{1,1}(nx:end,ny))-1, gh1{i}{1,1}(nx:end,ny),'color', [0.9 0.9 0.9], LineWidth=1)
    hold on
    gmean(:,i) = gh1{i}{1,1}(nx:end,ny);
end
plot(0:length(mean(gmean,2))-1,mean(gmean,2), '-k',LineWidth=1.5)
hold on
plot(0:length(quantile(gmean,0.05,2))-1,quantile(gmean,0.05,2),'--k',LineWidth=1)
hold on
plot(0:length(quantile(gmean,0.95,2))-1,quantile(gmean,0.95,2),'--k',LineWidth=1)
hold on
h1 = plot(0:length(gh{1,1}(nx:end,ny))-1,gh{1,1}(nx:end,ny),'-g',LineWidth=2);
hold on
h2 = plot(0:length(gh{1,1}(nx,ny:end))-1,gh{1,1}(nx,ny:end),'-r',LineWidth=2);
legend([h1 h2],{'Vertical', 'Horizontal'})
xlabel('Distance')
ylabel('Rank asymmetry value')
fontsize(gca, 14, 'points')
xlim([0 floor(nx/4)])
ylim([-0.0495 0.0495])
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')


%% Compute the variogram
icode = 1;
Z(Z==0)=0.0000001;
Z(Z==1)=0.9999999;

[gh, nh] = GeoStatFFT(x0, norminv(Z), icode);

% Simulation of random fields
[z]=fftma([4 50 50 50 0 0 0],1,45124241,100,nx,1,ny,1);
for i =1:100
    gh1{i} = GeoStatFFT(x0, z(:,i), icode);
end

figure(5)
for i=1:100
    plot(0:length(gh1{i}{1,1}(nx:end,ny))-1, gh1{i}{1,1}(nx:end,ny),'color', [0.9 0.9 0.9], LineWidth=1)
    hold on
    gmean(:,i) = gh1{i}{1,1}(nx:end,ny);
end
plot(0:length(mean(gmean,2))-1,mean(gmean,2), '-k',LineWidth=1.5)
hold on
plot(0:length(quantile(gmean,0.05,2))-1,quantile(gmean,0.05,2),'--k',LineWidth=1)
hold on
plot(0:length(quantile(gmean,0.95,2))-1,quantile(gmean,0.95,2),'--k',LineWidth=1)
hold on
h1 = plot(0:length(gh{1,1}(nx:end,ny))-1,gh{1,1}(nx:end,ny),'-g',LineWidth=2);
hold on
h2 = plot(0:length(gh{1,1}(nx,ny:end))-1,gh{1,1}(nx,ny:end),'-r',LineWidth=2);
legend([h1 h2],{'Vertical', 'Horizontal'})
xlabel('Distance')
ylabel('Variogram value')
fontsize(gca, 14, 'points')
xlim([0 floor(nx/4)])
ylim([0 1.4])
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
