# SpatialStatisticsFFT

This repository contains the computer codes associated with the manuscript: *Efficient Computation on Large Regular Grids of Higher-Order Spatial Statistics via Fast Fourier Transform*.

### MATLAB Live Script
We provide MATLAB Live Scripts that explain all the examples in detail and produce the figures displayed in the manuscript. The MATLAB Live Scripts need to be run in the same folder where the code is located. Codes are included in the example files. 
The following Live Scripts execute all the required code to generate the figures shown in the manuscript, along with additional outputs:
- `E1_GaussianCase.mlx` : 4.1. Synthetic Gaussian Fields
- `E2_Dunes.mlx` : 4.2. Dunes Dataset
- `E3A_AsymmetryFieldDir.mlx` : 4.3. Synthetic Non-Gaussian Fields
- `E3B_AsymmetryFieldRank.mlx` : 4.3. Synthetic Non-Gaussian Fields
- `E4_Chloride.mlx` : 4.4. Chloride concentrations in Baden-Wuerttemberg
- `E5_Xrays.mlx` : 4.5. X-ray absorption measurements of a soil column
- `E6_RhineRiver.mlx` : 4.6. One-dimensional example: Time-series analysis

### Example
The folder `Example` contains all the code related to the results presented in the manuscript, organized by number. Figures from the manuscript are stored in `.eps` format in their respective folders, and datasets are also included.

Gaussian random fields and synthetic datasets are generated using the FFTMA algorithm provided in this folder. The script `covar.m` is used to generate the covariance matrix, while `grille2.m` and `grille3.m` create 2D and 3D grids, respectively. The function `ECDF.m` computes the empirical cumulative distribution using ranking.

Additionally, `importfile.m` and `imagesc3d.m` are functions designed to import and plot the X-Ray dataset (Example 5).

The following scripts execute all the required code to generate the figures shown in the manuscript, along with additional outputs:
- `E1_GaussianCase.m` : 4.1. Synthetic Gaussian Fields
- `E2_Dunes.m` : 4.2. Dunes Dataset
- `E3A_AsymmetryFieldDir.m` : 4.3. Synthetic Non-Gaussian Fields
- `E3B_AsymmetryFieldRank.m` : 4.3. Synthetic Non-Gaussian Fields
- `E4_Chloride.m` : 4.4. Chloride concentrations in Baden-Wuerttemberg
- `E5_Xrays.m` : 4.5. X-ray absorption measurements of a soil column
- `E6_.RhineRiverm` : 4.6. One-dimensional example: Time-series analysis
  
The scripts `Time_GeoStatFFT.m` and `Time_Spatial_Spectral.m` perform time comparisons between several spatial statistics using the FFT method and a spatial program (`varioexp2d.m`), respectively.
  
## PROGRAMS

### `GeoStatFFT.m`
This function computes direct- and cross- spatial statistics in `nD` for up to `nvar` variables.

### `GeoStatFFT_ndir.m`
This function post-processes the output of `GeoStatFFT` to compute experimental directional or omnidirectional direct- and cross- spatial statistics.

#### Available Spatial Statistics:
- **1**: Variograms and cross-variograms
- **2**: Covariograms and cross-covariograms
- **3**: Variograms and pseudo-cross-variograms
- **4**: Centered covariances and cross-covariances (mean computed for the entire field)
- **5**: Non-centered covariances and cross-covariances (bivariate probabilities, for categorical data)
- **6**: Transiograms (for categorical data)
- **7**: Non-ergodic transiograms (for categorical data)
- **8**: Directional asymmetry and cross-directional asymmetry (Bárdossy & Hörning, 2017)
- **9**: Rank asymmetry and cross-rank asymmetry (Guthke, 2013)
- **10**: Rank correlation and cross-rank correlation (Bárdossy & Hörning, 2017)
- **11**: Third-order cumulant of a zero-mean random function (Dimitrakopoulos et al., 2010)

#### Modifications:
- Original version (`variofft2D.m`) written by **Denis Marcotte** (denis.marcotte@polymtl.ca), 1996. Developped spatial statistics **1**, **2**, **3**, **4**, **5**. 
- Adapted to `3D` by **Dimitri D'Or** (Ephesia Consult) on 2014/11/17, with additional capabilities for categorical data.  Added spatial statistics **6**, **7**. 
- Adapted by **Dany Lauzon** (Polytechnique Montréal, dany.lauzon@polymtl.ca) for `nD`, `nvar`, and non-collocated data. Added spatial statistics **8**, **9**, **10**, **11**. 

#### Reference:
- Marcotte, D., 1996. *Fast Variogram Computation with FFT*. Computers & Geosciences, 22(10), 1175-1186.
- Lauzon, D. and Hörning, S. *Efficient Computation on Large Regular Grids of High-Order Spatial Statistics via Fast Fourier Transform*. (In review)

---

