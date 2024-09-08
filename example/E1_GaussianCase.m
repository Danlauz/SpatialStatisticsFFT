%% Clear temporary variables
clear variables

%% Simulation 3 variables with 100 random fields

%Isotropic spherical model with range of 100m and unit sill (for covar).
model = [4 100 100 0 ]; c = 1;

% Generating a grid of size nx X ny
nx = 512;  
ny = 512;
x0 = grille2(1,nx,1,1,ny,1);

% Number of simulation
nbsim = 100;

% Simulating fields with FFTMA  
seed1 = 5451;
seed2 = 8964;
seed3 = 84;
[Z1]=fftma(model,c, seed1,nbsim, nx,1, ny,1);
[Z2]=fftma(model,c, seed2,nbsim, nx,1, ny,1);
[Z3]=fftma(model,c, seed3,nbsim, nx,1, ny,1);
 
% Theorical covariance function (for analysis and comparison)
[ktheo]=covar(1,(1:250)',model(1:2),1);
%% Linear Coregionalization model
% Build the simplify LCM
Y(:,1,:) = Z1 ;
Y(:,2,:) = 0.5 * Z1 + sqrt(1-0.5^2) * Z2;
Y(:,3,:) = Z3 ;
clear Z1 Z2 Z3
%% Plot the first LCM realizations
figure(1)
imagesc(reshape(Y(:,1,1),[nx ny]));
colorbar();
colormap 'jet'
clim([-3.2 3.2])
set(gca,'YDir','normal')
fontsize(gca, 12, 'points')

figure(2)
imagesc(reshape(Y(:,2,1),[nx ny]));
colorbar();
colormap 'jet'
clim([-3.2 3.2])
set(gca,'YDir','normal')
fontsize(gca, 12, 'points')

figure(3)
imagesc(reshape(Y(:,3,1),[nx ny]));
colorbar();
colormap 'jet'
clim([-3.2 3.2])
set(gca,'YDir','normal')
fontsize(gca, 12, 'points')
%% Compute the variograms and cross-variograms using FFT
icode = 1;
categ = 0; 
display = 1;

%Compute the nbsim variograms and cross-variograms of the LCM
for i=1:nbsim
[gh{i}, nh]=GeoStatFFT(x0,Y(:,:,i),icode);
end

% Plot the results
figure(4)
for k=1:3
    for j=1:3
        subplot(3,3,j+3*(k-1));
        means = 0;
        for i=1:100
            means = means + gh{i}{k,j}(nx:end,ny)/nbsim;
            plot(gh{i}{k,j}(nx:end,ny),'-','Color',[0.9 0.9 0.9])
            hold on
        end
        plot(means,'-k','LineWidth',2)
        if k==j
            hold on
            plot(1-ktheo,'-g','LineWidth',1.5)
        elseif (k==1 && j==2) || (j==1 && k==2)
            hold on
            plot(0.5*(1-ktheo),'-g','LineWidth',2)
        else
            hold on
            plot(0*(1-ktheo),'-g','LineWidth',2)
        end
        title(['Z',num2str(k),' vs. Z',num2str(j)]);
        xlim([0, round(2*nx/4)]);
    end
end

%% Compute the centered covariances and cross-covariances
categ = 0;
display = 1;
icode = 4;

%Compute the nbsim centered covariances and cross-covariances of the LCM
for i=1:nbsim
[gh{i}, nh]=GeoStatFFT(x0,Y(:,:,i),icode);
end

% Plot the results
figure(5)
for k=1:3
    for j=1:3
        subplot(3,3,j+3*(k-1));
        means = 0;
        for i=1:100
            means = means + gh{i}{k,j}(nx:end,ny)/nbsim;
            plot(gh{i}{k,j}(nx:end,ny),'-','Color',[0.9 0.9 0.9])
            hold on
        end
        plot(means,'-k','LineWidth',2)
        if k==j
            hold on
            plot(ktheo,'-g','LineWidth',1.5)
        elseif (k==1 && j==2) || (j==1 && k==2)
            hold on
            plot(0.5*(ktheo),'-g','LineWidth',2)
        else
            hold on
            plot(0*(ktheo),'-g','LineWidth',2)
        end
        title(['Z',num2str(k),' vs. Z',num2str(j)]);
        xlim([0, round(2*nx/4)]);
    end
end

%% Compute the non-centered covariances and cross-covariances
categ = 0;
display = 1;
icode = 5;

%Compute the nbsim non-centered covariances and cross-covariances of the LCM
for i=1:nbsim
[gh{i}, nh]=GeoStatFFT(x0,Y(:,:,i),icode);
end

% Plot the results
figure(6)
for k=1:3
    for j=1:3
        subplot(3,3,j+3*(k-1));
        means = 0;
        for i=1:100
            means = means + gh{i}{k,j}(nx:end,ny)/nbsim;
            plot(gh{i}{k,j}(nx:end,ny),'-','Color',[0.9 0.9 0.9])
            hold on
        end
        plot(means,'-k','LineWidth',2)
        if k==j
            hold on
            plot(ktheo,'-g','LineWidth',1.5)
        elseif (k==1 && j==2) || (j==1 && k==2)
            hold on
            plot(0.5*(ktheo),'-g','LineWidth',2)
        else
            hold on
            plot(0*(ktheo),'-g','LineWidth',2)
        end
        title(['Z',num2str(k),' vs. Z',num2str(j)]);
        xlim([0, round(2*nx/4)]);
    end
end

%% Compute the directional asymetry
categ = 0;
display = 1;
rank = 1; % We transform the data to the rank distribution (ECDF)
icode = 8;

%Compute the nbsim directional asymmetry of the LCM
for i=1:100
[gh{i}, nh]=GeoStatFFT(x0,Y(:,:,i),icode, 0, 0,rank);
end

figure(7)
for k=1:3
    for j=1:3
        subplot(3,3,j+3*(k-1));
        means = 0;
        for i=1:100
            means = means + gh{i}{k,j}(nx:end,ny)/100;
            plot(gh{i}{k,j}(nx:end,ny),'-','Color',[0.9 0.9 0.9])
            hold on
        end
        plot(means,'-k','LineWidth',2)
        hold on
        plot(0*(1-ktheo),'-g','LineWidth',2)

        title(['Z',num2str(k),' vs. Z',num2str(j)]);
        xlim([0, round(2*nx/4)]);
    end
end

%% Sampling 500 data (non-colocated) 

nbdata = 1000;
Z = Y(:,:,1); % Keep the first LCM realization

for i=1:3
    p = haltonset(2);
    id_x = floor(p(100+(nbdata*i +1):100+(nbdata*i +1)+nbdata-1,1)*nx)+1;
    id_y = floor(p(100+(nbdata*i +1):100+(nbdata*i +1)+nbdata-1,2)*nx)+1;
    id = sortrows([id_x id_y ],[1 2]);

    [ind{i}, I] = unique((id(:,1)-1)*(nx) + id(:,2));
    ind_all = (1:nx*ny)';


    figure(10+i)
    h = imagesc(reshape(Z(:,i),[nx ny]));
    hold on
    plot(id(:,1), id(:,2),'ok')
    colorbar();
    colormap 'jet'
    clim([-3.2 3.2])
    set(gca,'YDir','normal')
    xlim([0 nx])
    ylim([0 ny])

    Z(~ismember(ind_all,ind{i}),i)=nan;

    figure(100+i)
    scatter(id(I,1),id(I,2), 25, Z(ismember(ind_all,ind{i}),i), 'filled')
    colorbar();
    colormap 'jet'
    clim([-3.2 3.2])
    xlim([0 nx])
    ylim([0 ny])
end

[gh, nh] = GeoStatFFT(x0, Z, 5) ;

nbdist= 10;
ang=0;
tol_ang = 360;
max_dist = max(nx,ny)/3;
dist = [(0:nbdist-1);(1:nbdist)]'*(max_dist/nbdist);
[gh_omni, nh_omni, lag_omni] = GeoStatFFT_ndir(gh, nh, dist, ang, tol_ang);

nbdist= 10;
ang=[0,1,2,3,4,5,6,7]'*22.5;
tol_ang = 5;
max_dist = max(nx,ny)/3;
dist = [(0:nbdist-1);(1:nbdist)]'*(max_dist/nbdist);

[gh_dir, nh_dir, lag_dir] = GeoStatFFT_ndir(gh, nh, dist, ang, tol_ang);
%% 
dir = ang;
for i = 1:3
    for j = 1:3
        figure(40)
        subplot(3,3, ((i-1)*3 + j))
        plot(lag_omni{i,j}(:,1), gh_omni{i,j}(:,1), '*k')

        hold on
        if i == j
            plot(ktheo, '-g', 'LineWidth', 1.5)
        elseif (i == 1 && j == 2) || (j == 1 && i == 2)
            plot(0.5 * ktheo, '-g', 'LineWidth', 2)
        else
            plot(0 * ktheo, '-g', 'LineWidth', 2)
        end

        t = text(lag_omni{i,j}(:,1) + 0.05, gh_omni{i,j}(:,1) + 0.05, string(nh_omni{i,j}(:,1)));
        % Adjust font size for each text element
        for kk = 1:length(t)
            t(kk).FontSize = 10;
            t(kk).Rotation = 30;
        end
        xlim([0, round(max_dist)]);
        ylim([-0.3, 1.2]);
        title(['Covar. ', num2str(i), ' - ', num2str(j)]);
        xlabel('Distance')
        ylabel('Statistics value')
    end
end


figure(51)
i =3; j = 3;
for k=1:8
    subplot(3, 3, k)
    plot(lag_dir{i,j}(:,k), gh_dir{i,j}(:,k), '*k')

    hold on
    if i == j
        plot(ktheo, '-g', 'LineWidth', 1.5)
    elseif (i == 1 && j == 2) || (j == 1 && i == 2)
        plot(0.5 * ktheo, '-g', 'LineWidth', 2)
    else
        plot(0 * ktheo, '-g', 'LineWidth', 2)
    end

    t = text(lag_dir{i,j}(:,k) + 0.05, gh_dir{i,j}(:,k) + 0.05, string(nh_dir{i,j}(:,k)));

    % Adjust font size for each text element
    for kk = 1:length(t)
        t(kk).FontSize = 10;
        t(kk).Rotation = 60;
    end

    xlim([0, round(max_dist)]);
    ylim([-0.3, 1.2]);
    title(['Dir. ', num2str(dir(k))]);
    xlabel('Distance')
    ylabel('Statistics value')
end
