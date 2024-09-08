function [F] = ECDF(z)
% Function to transform data into univariate distribution values
% Syntax: [F] = ECDF(z)
% z is a matrix of size n x 1
%
% This function computes the empirical cumulative distribution function (ECDF)
% for the input matrix z and returns the transformed values in F. The function
% handles missing data (NaN) by assigning them a value of 0 in the output.
%
% Note: This function does not handle identical values.

%%% BE CAREFULL %%%
%%% WE break equalies randomly %%%

if any(isnan(z), 'all')
    F = z;
    idxnan = ~isnan(z);
    zvals = z(idxnan);
    zvals = zvals + rand(size(zvals))*10^-8;
    Fvals = zvals;

    [~, idx] = sort(zvals) ;
    rank =(0:length(zvals)-1)/(length(zvals)-1);
    Fvals(idx)=rank;
    F(idxnan) = Fvals; 
else
    F = z;
    z = z + rand(size(z))*10^-8;
    [~, idx] = sort(z) ;
    rank =(0:length(z)-1)/(length(z)-1);
    F(idx)=rank;
end


