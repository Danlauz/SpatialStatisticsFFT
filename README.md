# SpatialStatisticsFFT
Computer codes relate to the manuscript: Efficient Computation on Large Regular Grids of Higher-Order Spatial Statistics via Fast Fourier Transform

Example: The folder 'Example' contains all the code related to the results presented in the manuscript by number.

PROGRAMS
GeoStatFFT.m
% function [gh,nh]=GeoStatFFT(x0,z,icode,categ,display);
% Function to compute spatial statistics in nD for up to nvar variables.
% List of spatial statistics
%                      =1 : variograms and cross-variograms;
%                      =2 : covariograms and cross-covariograms
%                      =3 : variograms and pseudo-cross-variograms
%                      =4 : centered covariances and cross-covariances  (mean is computed for the whole field instead of according to the lags)
%                      =5 : non-centered covariances and cross-covariances (bivariate probabilities, for categorical data)
%                      =6 : transiograms (for categorical data)
%                      =7 : non-ergodic transiograms (for categorical data)
%                      =8 : Directional asymmetry and cross-directional asymmetry (Bardossy and Horning, 2017)
%                      =9 : Rank asymmetry and cross-rank asymmetry  (Guthke, 2013)
%                      =10: Rank correlation and cross-rank correlation (Bardossy and Horning, 2017)
%                      =11: Third-order cumulant of a zero mean random function  (Dimitrakopoulos et al., 2010)­
%
% Modification - Original version (1966)  variofft2D.m written by D. Marcotte, denis.marcotte@polymtl.ca
% Author: Dany Lauzon - Polytechnique Montréal
%                      Adapted to nD, nvar, and non-collocated data.
%                      =8 : Directional asymmetry and cross-directional asymmetry (Bardossy and Horning, 2017)
%                      =9 : Rank asymmetry and cross-rank asymmetry  (Guthke, 2013)
%                      =10: Rank correlation and cross-rank correlation (Bardossy and Horning, 2017)
%                      =11: Third-order cumulant of a zero mean random function  (Dimitrakopoulos et al., 2010)­
%  
% Author: Dimitri D'Or - Ephesia Consult - 2014/11/17 :
%                       Adapted to 3D.
%                       5 : bivariate probabilities, for categorical data
%                       6 : transiograms, for categorical data
%                       7 : non-ergodic transiograms, for categorical data
%
% Reference :
% Marcotte D., 1996. Fast Variogram computation with FFT. Computers & Geoscience, 22, 10, 1175-1186.
% Lauzon D. and Horning S. Efficient Computation on Large Regular Grids of High-Order Spatial Statistics via Fast Fourier Transform. (In review) 

GeoStatFFT_ndir
% function [gh_ndir, nh_ndir, lag_ndir] = GeoStatFFT_ndir(gh, nh, nbins, ndir);
% Function to post-process GeoStatFFT (Experimental directional or omnidirectional spatial statistics).
