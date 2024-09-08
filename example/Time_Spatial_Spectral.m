%% Variable number of data sampled
size = 128;
sample = [0.01 0.025 0.05 0.075 0.10 0.2 0.4 0.6 0.8 1];
nbsim = 10;

t1 = zeros(nbsim,length(sample));
t2 = zeros(nbsim,length(sample));
t3 = zeros(nbsim,length(sample));
for k = 1:length(size)
    for j = 1:length(sample)

        % Generate 100 fields
        nx = size(k);
        seed = 45152 ;
        model = [4 nx/4 nx/4 0 ];

        [Z] = fftma(model,1,seed,nbsim,nx,1,nx,1);
        x0 = grille2(1,nx,1,1,nx,1);
        % Sampling data
        nbdata = ceil((nx^2)*sample(j));
        sampling = true;
        if sampling
            p = haltonset(2);
            id_x = floor(p(100:100+nbdata*5,1)*nx)+1;
            id_y = floor(p(100:100+nbdata*5,2)*nx)+1;
            id = sortrows([id_x id_y ],[1 2]);

            [ind, I] = unique((id(:,1)-1)*(nx) + id(:,2));
            ind_all = (1:nx*nx)';
            ind= ind(1:nbdata);
        end
        % 1 simulations
        for i=1:nbsim
            [i,j,k]
            z = Z(:,i);

            % Spatial, All points
            x = [x0, z];
            x = x(ismember(ind_all,ind),:);
            nbclas=nx;
            lclas=[(0:nbclas-1);(1:nbclas)]'*(nx/nbclas);
            vdir=[0]';
            vreg=ones(1,1)*180;
            tic
            [gexp]=varioexp2d(x,nbclas,lclas,vdir,vreg);
            t1(i,j) = toc ;

            % Spatial, 8 lags, max nx/2
            x = [x0, z];
            x = x(ismember(ind_all,ind),:);
            nbclas=10;
            lclas=[(0:nbclas-1);(1:nbclas)]'*(nx/nbclas);
            vdir=[0,1,2,3,4,5,6,7]'*22.5;
            vreg=ones(8,1)*5;

            tic
            [gexp]=varioexp2d(x,nbclas,lclas,vdir,vreg);
            t2(i,j) = toc ;

            % FFT, all points, processing to 8 lags, max nx/2
            nbins = 10;
            ndir = 16; % While the variogram is symmetric, this does not hold for spatial asymmetry.
            tic
            [gh, nh]=GeoStatFFT(x0, z, 1);
            nbdist= 20;
            ang=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]'*22.5;
            tol_ang = 10;
            max_dist = nx/2;
            dist = [(0:nbdist-1);(1:nbdist)]'*(max_dist/nbdist);
            [gh_omni, nh_omni, lag_omni] = GeoStatFFT_ndir(gh, nh, dist, ang, tol_ang);
            t3(i,j) = toc ;

            [t1(i,j) t2(i,j) t3(i,j)]
        end
    end
end
t1 = mean(t1,1);
t2 = mean(t2,1);
t3 = mean(t3,1);

figure(1)
loglog( sample*size.^2, t1, 'og','MarkerEdgeColor','green', 'MarkerFaceColor','green')
hold on
loglog( sample*size.^2, t2, 'or','MarkerEdgeColor','red', 'MarkerFaceColor','red' )
hold on
loglog( sample*size.^2, t3, 'ob','MarkerEdgeColor','blue', 'MarkerFaceColor','blue' )
grid
fontsize(gca, 12,'points')
xlim([100 20000])
ylim([0.0005 50])
xlabel('Number of data sampled ','FontSize',16)
ylabel('Time (s)','FontSize',16)
legend('Variogram (All Points)', 'Variogram (8 Directions)', 'Variogram (FFT)', 'location', 'northwest')

%% Variable number of grid points with N_x = N_y
size = [2 4 8 16 32 64 128 256 512 1024 2048];
nbsim = 10;

t1 = zeros(nbsim,length(size));
t2 = zeros(nbsim,length(size));
t3 = zeros(nbsim,length(size));
for k = 1:length(size)

    % Generate 100 fields
    nx = size(k);
    seed = 45152 ;
    model = [4 nx/4 nx/4 0 ];

    [Z] = fftma(model,1,seed,nbsim,nx,1,nx,1);
    x0 = grille2(1,nx,1,1,nx,1);
    % Sampling data
    nbdata = ceil((nx^2));
    sampling = true;
    if sampling
        p = haltonset(2);
        id_x = floor(p(100:100+nbdata*5,1)*nx)+1;
        id_y = floor(p(100:100+nbdata*5,2)*nx)+1;
        id = sortrows([id_x id_y ],[1 2]);

        [ind, I] = unique((id(:,1)-1)*(nx) + id(:,2));
        ind_all = (1:nx*nx)';
        ind= ind(1:nbdata);
    end
    % 1 simulations
    for i=1:nbsim
        [i,k]
        z = Z(:,i);
        if k<8
            % Spatial, All points
            x = [x0, z];
            x = x(ismember(ind_all,ind),:);
            nbclas=nx;
            lclas=[(0:nbclas-1);(1:nbclas)]'*(nx/nbclas);
            vdir=[0]';
            vreg=ones(1,1)*180;
            tic
            [gexp]=varioexp2d(x,nbclas,lclas,vdir,vreg);
            t1(i,k) = toc ;
        end

        if k<8
            % Spatial, 8 lags, max nx/2
            x = [x0, z];
            x = x(ismember(ind_all,ind),:);
            nbclas=10;
            lclas=[(0:nbclas-1);(1:nbclas)]'*(nx/nbclas);
            vdir=[0,1,2,3,4,5,6,7]'*22.5;
            vreg=ones(8,1)*5;

            tic
            [gexp]=varioexp2d(x,nbclas,lclas,vdir,vreg);
            t2(i,k) = toc ;

        end

        % FFT, all points, processing to 8 lags, max nx/2
        nbins = 10;
        ndir = 16; % While the variogram is symmetric, this does not hold for spatial asymmetry.
        tic
        [gh, nh]=GeoStatFFT(x0, z, 1);
        nbdist= 20;
        ang=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]'*22.5;
        tol_ang = 10;
        max_dist = nx/2;
        dist = [(0:nbdist-1);(1:nbdist)]'*(max_dist/nbdist);
        [gh_omni, nh_omni, lag_omni] = GeoStatFFT_ndir(gh, nh, dist, ang, tol_ang);
        t3(i,k) = toc ;

        [t1(i,k) t2(i,k) t3(i,k)]
    end
end
t1 = mean(t1,1);
t2 = mean(t2,1);
t3 = mean(t3,1);


figure(2)
loglog( size.^2, t1, 'og','MarkerEdgeColor','green', 'MarkerFaceColor','green')
hold on
loglog( size.^2, t2, 'or','MarkerEdgeColor','red', 'MarkerFaceColor','red' )
hold on
loglog( size.^2, t3, 'ob','MarkerEdgeColor','blue', 'MarkerFaceColor','blue' )
grid
fontsize(gca, 12,'points')
xlim([1 10000000])
ylim([0.0001 100])
xlabel('Number of grid points with N_x = N_y','FontSize',16)
ylabel('Time (s)','FontSize',16)
legend('Variogram (All Points)', 'Variogram (8 Directions)', 'Variogram (FFT)', 'location', 'northwest')
