function [gh,nh]=GeoStatFFT(x0,z,icode,categ,display,rank)
% function [gh,nh]=GeoStatFFT(x0,z,icode,categ,display);
% Function to compute spatial statistics in nD for up to nvar variables.
%
% Missing data are replace by zeros. A data are on a (possibly incomplete) regular grid.
% The program computes spatial statistics in the frequency domain using ND-FFT.
%
% INPUT :
% x0       n by d     matrix of coordinates (regular grid)
%
% z        n by nvar  matrix of values for the variables. Each line is
%                     associated with the corresponding line vector of
%                     coordinates in the x0 matrix, and each column corresponds
%                     to a different variable.
%                     Missing values are indicated by NaN. We replace NaN values with
%                     zeros in the code.
%
% icode               a code to indicate which function to compute
%                      =1 : variograms and cross-variograms;
%                      =2 : covariograms and cross-covariograms
%                      =3 : variograms and pseudo-cross-variograms
%                      =4 : centered covariances and cross-covariances  (mean is computed for the whole field instead of according to the lags)
%                      =5 : non-centered covariances and cross-covariances (bivariate probabilities, for categorical data)
%                      =6 : transiograms (for categorical data)
%                      =7 : non-ergodic transiograms (for categorical data)
%                      =8 : Directional asymmetry (Bardossy and Horning, 2017)
%                      =9 : Bivariate asymmetry (Guthke, 2013)
%                      =10: Rank correlation and cross-rank correlation 
%                      =11: Third-order cumulant of a zero mean random function  (Dimitrakopoulos et al., 2010)­
%
% categ    boolean    tells if the variables are categorical (set to 1) or not
%                     (set to 0, default)
% display  boolean    tells if a plot must be displayed (set to 1) or not
%                     (set to 0, default)
% rank     boolean    tells if we need to perform rank transformation (set to 1) or not
%                     (set to 0, default)

% OUTPUT :
% gh       nvar by nvar cell array of nx by ny (direct- and cross-) maps for
%                     variables i and j depending on icode.
% nh       nnvar by nvar cell array of nx by ny (direct- and cross-) maps of number of pairs
%                     available to compute the structural function.
%
% This program uses the functions FFTn, IFFTn, FFTSHIFT and CONJ which are standard MATLAB functions.
%
% Modification - Original version  variofft2D.m written by D. Marcotte, denis.marcotte@polymtl.ca
% Author: Dany Lauzon - Polytechnique Montréal
%                      =8 : Directional asymmetry (Bardossy and Horning, 2017)
%                      =9 : Bivariate asymmetry (Guthke, 2013)
%                      =10: Rank correlation and cross-rank correlation 
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

if nargin<6
    rank=0;
end
if nargin<5
    display=0;
end
if nargin<4
    categ=0;
end


% Number of variables (or fields)
[n_nodes,nvar]=size(z);
if categ
    nvar=length(unique(z(~isnan(z))));
end

%%% Finding the parameters of the grid
minc=min(x0);
maxc=max(x0);

dc = zeros(size(x0,2),1);
for i=1:size(x0,2)
    dum1=diff(unique(x0(:,i)));
    if isempty(dum1)
        dc(i)=1;
    else
        dc(i)=dum1(1);
    end
end
nc=((maxc-minc)./dc')+1;

%%% Reformatting the data if categorical
if categ
    idnan=isnan(z);
    zc=zeros(n_nodes,nvar);
    for i=1:nvar
        zc(z==i,i)=1;
        zc(idnan,i)=nan;       % to keep the nan in the computation
    end
    z=zc;
    clear zc;
end

%%% Initialization
Z=cell(nvar,1);
Zid=cell(nvar,1);

if length(nc)==1
    nc=[nc 1];
end
for i=1:nvar
    Z{i}=reshape(z(:,i),nc);
end

gh=cell(nvar,nvar);
[n,p,q]=size(Z{1});          % dimensions of data matrix

% find the closest multiple of 8 to obtain a good compromise between
% speed (a power of 2) and memory required

nrows=2*n-1;
ncols=2*p-1;
nz=2*q-1;
nr2=ceil(nrows/2)*2;
nc2=ceil(ncols/2)*2;
nz2=ceil(nz/2)*2;
nv=[nr2,nc2,nz2];


% form an indicator  matrix:  1's for all data values
%                             0's for missing values
% in data matrix, replace missing values by 0;

for i=1:nvar
    Zid{i}=~isnan(Z{i});                        % 1 for a data value; 0 for missing
    Z{i}(~Zid{i})=0;                            % missing replaced by 0
end

% Apply the rank transformation
Fz= cell(nvar,1);
if rank ==1
    for i=1:nvar
        Fz{i}=ECDF(z(:,i));
        Fz{i} = reshape(Fz{i},nc);
        Fz{i}(~Zid{i})=0; 
    end
else
    for i=1:nvar
        Fz{i}=z(:,i);
        Fz{i} = reshape(Fz{i},nc);
        Fz{i}(~Zid{i})=0; 
    end
end

% Preparation
fx=cell(nvar,1);
fx2=cell(nvar,nvar);
fxid=cell(nvar,nvar);

% Compute probability of each facies (if categ)
if icode(1)==6 && ~categ
    error('Transiograms are for categorical data, categ = 1 not 0')
end
if icode(1)==6 && categ
    prop = zeros(nvar,1);
    for i=1:nvar
        prop(i)=nansum(z(:,i))/(length(z(:,i))-sum(idnan));
    end
end

% Compute the mean if data are centered
if icode(1)==4 || icode(1)==11 
    m = zeros(nvar,1);
    for i=1:nvar
        m(i)=sum(sum(Z{i}(Zid{i})))/sum(sum(Zid{i}));
        Z{i}(Zid{i})=Z{i}(Zid{i})-m(i);
    end
end

% Compute the Fourier transform of variables and indicators
for i=1:nvar
    for j=1:nvar
        if i==j
            fx{i}=fftn(Z{i},nv);                         % fourier transform of Z{i}
            fxid{i,i}=fftn(Zid{i},nv);
        else
            fxid{i,j}=fftn(Zid{i}.*Zid{j},nv);
        end
        fx2{i,j}=fftn(Z{i}.*Z{j},nv);                    % fourier transform of Z{i}*Z{j}
    end
end

% Compute number of pairs
for i=1:nvar
    for j=1:nvar
        if icode(1)==1
            nh{i,j}=round(real(ifftn(conj(fxid{i,j}).*fxid{i,j})));  % number of pairs f(x), g(x), f(x+h), g(x+h)
        else
            nh{i,j}=round(real(ifftn(conj(fxid{i,i}).*fxid{j,j})));  % number of pairs f(x), g(x+h)
        end
    end
end

% compute the different structural functions according to icode
switch icode(1)
    case 1 % variograms and cross-variograms are computed

        for i=1:nvar
            for j=1:nvar
                t1=fftn(Z{i}.*Zid{j},nv);
                t2=fftn(Z{j}.*Zid{i},nv);
                t12=fftn(Z{i}.*Z{j},nv);
                gh{i,j}=real(ifftn(conj(fxid{i,j}).*t12+conj(t12).*fxid{i,j}-conj(t1).*t2-t1.*conj(t2)))./max(nh{i,j},1)/2;
            end
        end

    case 2 % covariograms and cross-covariograms are computed

        for i=1:nvar
            for j=1:nvar
                m_tail=real(ifftn(conj(fx{i}).*fxid{j,j}))./max(nh{i,j},1); % computes the tail means
                m_head=real(ifftn(conj(fxid{i,i}).*fx{j}))./max(nh{i,j},1); % computes the head means

                gh{i,j}=real((ifftn(conj(fx{i}).*fx{j}))./max(nh{i,j},1)-m_tail.*m_head);
            end
        end

    case 3 % variograms and pseudo-cross-variograms are computed

        for i=1:nvar
            for j=1:nvar
                gh{i,j}=real(ifftn(fxid{j,j}.*conj(fx2{i,i})+conj(fxid{i,i}).*fx2{j,j}-2*conj(fx{i}).*fx{j}))./max(nh{i,j},1)/2;
            end
        end

    case {4,5}

        for i=1:nvar
            for j=1:nvar
                gh{i,j}=real((ifftn(conj(fx{i}).*fx{j}))./max(nh{i,j},1));
            end
        end

    case 6 % Transiograms
        for i=1:nvar
            for j=1:nvar
                gh{i,j}=(real((ifftn(conj(fx{i}).*fx{j}))./max(nh{i,j},1)))/prop(i);
            end
        end

    case 7  % non ergodic transiograms

        fx_all=0;
        for i=1:nvar
            fx_all=fx_all+fx{i};
        end

        for i=1:nvar
            propi=round(real(ifftn(conj(fx{i}).*fx_all)));
            for j=1:nvar
                gh{i,j}=real((ifftn(conj(fx{i}).*fx{j})))./max(propi,1);
            end
        end

    case 8  % Directional asymmetry (Bardossy and Horning, 2017)
        for i=1:nvar
            for j=1:nvar
                f3=fftn(Fz{i}.*Fz{i}.*Fz{i},nv);
                g3=fftn(Fz{j}.*Fz{j}.*Fz{j},nv);
                f2=fftn(Fz{i}.*Fz{i},nv);
                g2=fftn(Fz{j}.*Fz{j},nv);
                f1=fftn(Fz{i},nv);
                g1=fftn(Fz{j},nv);
                gh{i,j}=real(ifftn(conj(f3).*fxid{j,j}-3*conj(f2).*g1+3*conj(f1).*g2-conj(fxid{i,i}).*g3))./max(nh{i,j},1);
            end
        end

    case 9  % Bivariate asymmetry (Guthke, 2013)
        for i=1:nvar
            for j=1:nvar
                f3=fftn(Fz{i}.*Fz{i}.*Fz{i},nv);
                g3=fftn(Fz{j}.*Fz{j}.*Fz{j},nv);
                f2=fftn(Fz{i}.*Fz{i},nv);
                g2=fftn(Fz{j}.*Fz{j},nv);
                f1=fftn(Fz{i},nv);
                g1=fftn(Fz{j},nv);
                gh{i,j}=real(ifftn(      conj(f3).*fxid{j,j} + 3.*conj(f2).*g1        - 3.*conj(f2).*fxid{j,j} +...
                                      3.*conj(f1).*g2        - 6.*conj(f1).*g1        + 3.*conj(f1).*fxid{j,j} +...
                                         conj(fxid{i,i}).*g3 - 3.*conj(fxid{i,i}).*g2 + 3.*conj(fxid{i,i}).*g1 -conj(fxid{i,i}).*fxid{j,j}    ))./max(nh{i,j},1);
            end
        end

    case 10  % Rank correlation
        for i=1:nvar
            for j=1:nvar           
                f1=fftn(Fz{i},nv);
                f2=fftn(Fz{j},nv);
                gh{i,j}=12*real(ifftn( conj(f1).*f2 - 0.5*conj(f1).*fxid{j,j} - 0.5*conj(fxid{i,i}).*f2 + 0.25*conj(fxid{i,i}).*fxid{j,j} ))./max(nh{i,j},1);
            end
        end

    case 11  % High-order statistics, Cumulant tri-point (order 3)  (Article : Dimitrakopoulos (2009))
        clear nh
        h2i=icode(2); h2j=icode(3); %h2   fix to one direction

        Fzg=cell(1,nvar);
        Fzf=cell(1,nvar);
        for i=1:nvar
            %finding fonction f(x).*f(x+h2)=g(x,x+h2) and pairs existence
            id=zeros(3*(n-1)+1,3*(p-1)+1);
            id(n:1:2*n-1,p:1:2*p-1)=Z{i};
            g=id(n:1:2*n-1,p:1:2*p-1).*id(n+h2i:1:2*n-1+h2i,p+h2j:1:2*p-1+h2j);
            idh2=~(g==0);

            %Transform to ecdf
            fg{i}=g;
            ff{i}=Z{i};
            fg{i}(~idh2)=0;
            ff{i}(~Zid{i})=0;
            % finding number of tri-pairs (x,x+h1,x+h2)
            nh{i,i}=round(real(ifftn(conj(fftn(idh2,nv)).*fxid{i,i})));
            
            %  f(x)*f(x+h1)*f(x+h2)
            gh{i,i}=real((ifftn(conj(fftn(ff{i},nv)).*fftn(fg{i},nv))))./max(nh{i,i},1);
            
        end
end

% reduce matrices to required size
t=nv/2+1;
for i=1:nvar
    if icode(1) ~= 11
        seq = 1:nvar;
    else
        seq = i;
    end
    for j=seq

        if isempty(gh{i,j})
            continue;
        end
        ghtemp=fftshift(gh{i,j});
        gh{i,j}=ghtemp(t(1)-n+1:t(1)+n-1,t(2)-p+1:t(2)+p-1,t(3)-q+1:t(3)+q-1);
        nhtemp=fftshift(nh{i,j});
        nh{i,j}=nhtemp(t(1)-n+1:t(1)+n-1,t(2)-p+1:t(2)+p-1,t(3)-q+1:t(3)+q-1);
    end

end

% % Display graph
if display && (length(size(gh{1,1}))<3) % if display and we are in 2D
    figure
    for i=1:nvar
        for j=1:nvar
            if isempty(gh{i,j})
                continue;
            end
            subplot(nvar,nvar,j+nvar*(i-1));
            %      pcolor(Xmesh,Ymesh,gh{i,j});
            imagesc(gh{i,j}')
            hold on
            plot([nc(1) nc(1)],[0 nc(2)*2],'-k')
            hold on
            plot([0 nc(1)*2],[nc(2) nc(2)],'-k')

            if icode(1) == 1
                title(['Var. ',num2str(i),' - ',num2str(j)]);
            elseif icode(1) == 2
                title(['Covario. ',num2str(i),' - ',num2str(j)]);
            elseif icode(1) == 3
                title(['Pseudo-Var. ',num2str(i),' - ',num2str(j)]);
            elseif icode(1) == 4
                title(['Centered Covar. ',num2str(i),' - ',num2str(j)]);
            elseif icode(1) == 5 && categ == 0
                title(['Covar. ',num2str(i),' - ',num2str(j)]);
            elseif icode(1) == 5 && categ == 1
                title(['Bivar. Prob. ',num2str(i),' - ',num2str(j)]);
            elseif icode(1) == 6
                title(['Transio. ',num2str(i),' - ',num2str(j)]);
            elseif icode(1) == 7
                title(['Var. ',num2str(i),' - ',num2str(j)]);
            elseif icode(1) == 8
                title(['Dir. Asy. ',num2str(i),' - ',num2str(j)]);
            elseif icode(1) == 9
                title(['Biv. Asy. ',num2str(i),' - ',num2str(j)]);
            elseif icode(1) == 10
                title(['Rank Cor. ',num2str(i),' - ',num2str(j)]);
            elseif icode(1) == 11
                title(['Third-order cumulant ',num2str(i)]);
            end

            axis equal
            
            % Set the axis limits to center the image at (0, 0)
            xlim([round(3*nc(1)/4), round(5*nc(1)/4)]);
            ylim([round(3*nc(2)/4), round(5*nc(2)/4)]);

            shading flat;
            set(gca,'YDir','normal')

            if icode(1)==8 || icode(1)==9 || icode(1)==10 
                clim([-0.05 0.05]);
            elseif  categ == 0 && icode(1)==6
                clim([0 1]);
            elseif icode(1)==11
                clim([-0.05 0.05])
            elseif categ == 1 && icode(1)==6
                clim([0 0.75])
            elseif categ == 1 && icode(1)==5
                clim([0 0.4])
            else
                clim([0 1])
            end
            colorbar
        end
    end
end

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