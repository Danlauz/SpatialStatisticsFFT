%% Clear temporary variables
clear variables
nx = 250;
ny = 250;
x0 = grille2(1,nx,1,1,ny,1);

load ("E3 - SyntheticAsymmetry/Synthetic_assymetric_field_c0.mat")

Z = [Z_Ass, Z_familly(:,51)];

%%
% Rank asymmetry field
figure(1)
imagesc(reshape(Z(:,1),[nx ny]));
colorbar();
clim([-3.25 3.25])
colormap 'jet'
set(gca,'YDir','normal')

% One realization of the familly
figure(2)
imagesc(reshape(Z(:,2),[nx ny]))
colorbar();
clim([-3.25 3.25])
colormap 'jet'
set(gca,'YDir','normal')

% Ensemble of the familly use to generate the rank assymetry field
figure(3)
a=0;
for i = round(1:  (size(Z_familly,2)-1)/9 : size(Z_familly,2)-1)
    a=a+1;
    subplot(3,3,a)
    imagesc(reshape(Z_familly(:,i),[nx ny]))
    colorbar();
    clim([-3.25 3.25])
    colormap 'jet'
    axis equal
    xlim([0 nx]);
    ylim([0 ny]);
    set(gca,'YDir','normal')
end

% Multiple Realization of a set of parameter inside the familly
figure(4)
a=0;
for i = round(1:  (size(Z_familly,2)-1)/9 : size(Z_familly,2)-1)
    a=a+1;
    subplot(3,3,a)
    imagesc(reshape(Y(:,i),[nx ny]))
    colorbar();
    clim([-3.25 3.25])
    colormap 'jet'
    axis equal
    xlim([0 nx]);
    ylim([0 ny]);
    set(gca,'YDir','normal')
end
%% Compute the variograms and cross-variograms
categ = 0;
display = 1;
icode = 1;
[gh, nh]=GeoStatFFT(x0,Z,icode,categ,display);

figure(8)
% Simulation of random fields
[z]=fftma([1 0 0 0 0 0 0 ; 2 50/3 50/3 50/3 0 0 0],[0.13 ; 0.87],45124241,100,nx,1,ny,1);
for i =1:100
    gh1{i} = GeoStatFFT(x0,z(:,i),icode,categ,0);
end

figure(9)
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
legend([h1 h2],{'Vertical', 'Horizontal'},'FontSize',14)
fontsize(gca, 14, 'points')
xlabel('Distance','FontSize',16)
ylabel('Variogram value','FontSize',16)
xlim([0 floor(nx/4)])
ylim([0 1.4])
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
%% Compute the directional asymetry
categ = 0;
display = 1;
icode = 8;
rank = 1;
[gh, ~]=GeoStatFFT(x0,Z,icode,categ,display,rank);

figure(10)
% Simulation of random fields
[z]=fftma([6 50 50 50 0 0 0],1,45124241,100,nx,1,ny,1);
for i =1:100
    gh1{i} = GeoStatFFT(x0,z(:,i),icode,categ,0, rank);
end

figure(11)
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
legend([h1 h2],{'Vertical', 'Horizontal'},'FontSize',14)
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
rank = 1;
[gh, ~]=GeoStatFFT(x0,Z,icode,categ,display, rank);

figure(12)
% Simulation of random fields
[z]=fftma([6 50 50 50 0 0 0],1,45124241,100,nx,1,ny,1);
for i =1:100
    gh1{i} = GeoStatFFT(x0,z(:,i),icode,categ,0, rank);
end

figure(13)
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
legend([h1 h2],{'Vertical', 'Horizontal'},'FontSize',14)
xlabel('Distance')
ylabel('Rank asymmetry value')
fontsize(gca, 14, 'points')
xlim([0 floor(nx/4)])
ylim([-0.0495 0.0495])
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')