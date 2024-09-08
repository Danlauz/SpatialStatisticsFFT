# SpatialStatisticsFFT

This repository contains the computer codes associated with the manuscript: *Efficient Computation on Large Regular Grids of Higher-Order Spatial Statistics via Fast Fourier Transform*.

### Example
The folder `Example` contains all the code related to the results presented in the manuscript, organized by result number.

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

