function [gh_ndir, nh_ndir, lag_ndir] = GeoStatFFT_ndir(gh, nh, dist, ang, tol_ang)
% function [gh_ndir, nh_ndir, lag_ndir] = GeoStatFFT_ndir(gh, nh, nbins, ndir);
% Function to post-process GeoStatFFT (Experimental directional or omnidirectional spatial statistics).
%
% Missing data are replaced by zeros. The data are assumed to be on a (possibly incomplete) regular grid.
% The program computes spatial statistics in the frequency domain using ND-FFT.
%
% INPUT :
% gh       nvar by nvar cell array of nx by ny (direct- and cross-) maps for
%          variables i and j depending on icode.
% nh       nvar by nvar cell array of nx by ny (direct- and cross-) maps of the number of pairs
%          available to compute the structural function.
% nbdist  Number of distance class to compute the experimental spatial statistics
% max_dist Maximum distance for lags          
% ang      Number of angle to compute the experimental spatial statistics
% tol_ang  Tolerance on the angle.
%
% OUTPUT :
% gh_ndir  nvar by nvar cell array of nx by ny (direct- and cross-) maps for
%          variables i and j depending on icode.
% nh_ndir  nvar by nvar cell array of nx by ny (direct- and cross-) maps of the number of pairs
%          available to compute the structural function.
% lag_ndir nvar by nvar cell array of nx by ny (direct- and cross-) maps of distances
%          to compute the structural function.

% Number of variables
nvar = size(gh,1);

% Dimension of field
nc = (size(gh{1,1})-1)/2;
ndim = length(nc);

% Generate matrix of lags

if ndim == 2
    [X,Y] = meshgrid(-nc(2):nc(2),-nc(1):nc(1));
    lags = sqrt(X.^2 + Y.^2);
    angle = flip(rad2deg(atan(Y./X))) ;
    angle = angle + 180*(angle>0).*(X<0);
    angle = angle + 180*(angle<=0).*(Y<=0).*(X<0);
    angle = angle + 360*(angle<0).*(X>=0);
    angle(nc(1)+1, nc(2)+1)=0;
else
    error('ndim not equal 2')
end

tol_dist = dist;
tol_dist = [-1 0 ; tol_dist]; % artifical class to keep the result at lags 0.
tol_ang = [ang-tol_ang ang+tol_ang];

% Iteration over the number of variables
for i = 1 : nvar
    for j= 1 : nvar
        for k = 1 : size(tol_dist,1)
            for kk = 1 : size(tol_ang,1)
                if k==1
                    kkk=1;
                else
                    kkk=kk;
                end
                id =  (lags > tol_dist(k,1)) & (tol_dist(k,2) >= lags) & (angle > tol_ang(kkk,1)) & (tol_ang(kkk,2) >= angle) ;
                nh_ndir{i,j}(k,kk)  = sum( nh{i,j}(id));
                gh_ndir{i,j}(k,kk)  = sum((gh{i,j}(id) .* nh{i,j}(id)))/ nh_ndir{i,j}(k,kk);
                lag_ndir{i,j}(k,kk) = sum((lags(id) .* nh{i,j}(id)))/ nh_ndir{i,j}(k,kk);
            end
        end
    end
end
