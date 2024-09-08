size = [4 8 16 32 64 128 256 512 1024 2048];
sample = [1];

nbsim = 10;
size_dim = length(size);
t1 = zeros(nbsim ,size_dim ); t2 = zeros(nbsim ,size_dim ); t3 = zeros(nbsim ,size_dim );
t4 = zeros(nbsim ,size_dim ); t5 = zeros(nbsim ,size_dim ); t6 = zeros(nbsim ,size_dim );


for k = 1:length(size)
    for j = 1:length(sample)
        % Parameters for experimental spatial statitics
        nbdist= 20;
        ang=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]'*22.5;
        tol_ang = 10;
        max_dist = size(k)/2;
        dist = [(0:nbdist-1);(1:nbdist)]'*(max_dist/nbdist);

        % Generate nbsim  fields
        nx = size(k);
        seed = 45152 ;
        model = [4 nx/4 nx/4 0 ];

        [Z] = fftma(model,1,seed,nbsim ,nx,1,nx,1);
        x0 = grille2(1,nx,1,1,nx,1);
        % Sampling data
        nbdata = (nx^2)*sample;
        sampling = true;
        if sampling
            p = haltonset(2);
            id_x = floor(p(100 :100+nbdata*5,1)*nx)+1;
            id_y = floor(p(100:100+nbdata*5,2)*nx)+1;
            id = sortrows([id_x id_y ],[1 2]);

            [ind, I] = unique((id(:,1)-1)*(nx) + id(:,2));
            ind_all = (1:nx*nx)';
            ind= ind(1:nbdata);
        end
        % nbsim  simulations
        for i=1:nbsim 
            [i,j,k]
            z = Z(:,i);


            % FFT, all points, processing to 8 lags, max nx/2
            tic
            [gh, nh]=GeoStatFFT(x0, z, 1);
            [gh_omni, nh_omni, lag_omni] = GeoStatFFT_ndir(gh, nh, dist, ang, tol_ang);
            t1(i,k) = toc ;
            % FFT, all points, processing to 8 lags, max nx/2
            tic
            [gh, nh]=GeoStatFFT(x0, z, 4);
            [gh_omni, nh_omni, lag_omni] = GeoStatFFT_ndir(gh, nh, dist, ang, tol_ang);
            t2(i,k) = toc ;

            % FFT, all points, processing to 8 lags, max nx/2
            tic
            [gh, nh]=GeoStatFFT(x0, z, 5);
            [gh_omni, nh_omni, lag_omni] = GeoStatFFT_ndir(gh, nh, dist, ang, tol_ang);
            t3(i,k) = toc ;

            % FFT, all points, processing to 8 lags, max nx/2
            tic
            [gh, nh]=GeoStatFFT(x0, z, 8);
            [gh_omni, nh_omni, lag_omni] = GeoStatFFT_ndir(gh, nh, dist, ang, tol_ang);
            t4(i,k) = toc ;

            % FFT, all points, processing to 8 lags, max nx/2
            tic
            [gh, nh]=GeoStatFFT(x0, z, 9);
            [gh_omni, nh_omni, lag_omni] = GeoStatFFT_ndir(gh, nh, dist, ang, tol_ang);
            t5(i,k) = toc ;

            % FFT, all points, processing to 8 lags, max nx/2
            tic
            [gh, nh]=GeoStatFFT(x0, z, 10);
            [gh_omni, nh_omni, lag_omni] = GeoStatFFT_ndir(gh, nh, dist, ang, tol_ang);
            t6(i,k) = toc ;

            [t1(i,k) t2(i,k)]
        end
    end
end

%%

figure(1)
loglog( size.^2, mean(t1), 'or','MarkerEdgeColor','red', 'MarkerFaceColor','red')
hold on
loglog( size.^2, mean(t2), 'og','MarkerEdgeColor','green', 'MarkerFaceColor','green' )
hold on
loglog( size.^2, mean(t3), 'ob','MarkerEdgeColor','blue', 'MarkerFaceColor','blue' )
hold on
loglog( size.^2, mean(t4), 'oc','MarkerEdgeColor','cyan', 'MarkerFaceColor','cyan' )
hold on
loglog( size.^2, mean(t5), 'om','MarkerEdgeColor','magenta', 'MarkerFaceColor','magenta' )
hold on
loglog( size.^2, mean(t6), 'oy','MarkerEdgeColor','yellow', 'MarkerFaceColor','yellow' )
grid
fontsize(gca, 12,'points')
xlim([50 10000000])
ylim([0.001 30])
xlabel('Number of grid points with N_x = N_y','FontSize',16)
ylabel('Time (s)','FontSize',16)
legend('Variogram', 'Cen. Covariance','Non-Cen. Covariance', 'Dir. Asymmetry', 'Spa. Asymmetry', 'Rank Corr.', 'location','northwest')

