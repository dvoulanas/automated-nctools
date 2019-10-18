# automated-nctools

automated-nctools is a suite of automated MATLAB functions. The first on is the one below.

1. function [ncdata,selpoints,distselind]=ncextract(filename,ctlpoints,DistClass,k)

This function searches for the nearest k-neighbors to the control points, within XY or does
exact comparison to the grid cells of the gridded data coords of regular lonlat rectagular grid of netcdfs files or makes exact comparisons or each corresponding grid cell
filename is a cell vector with fullfile name of the netcdf fiel to be accessed
ctlpoints are Num_ctl_points by 2 matrix. Each row it's a pair of XY coordinates.
seldist is equal to the rows of ctlpoints, and the number of columns is k, which corresponds to the number of neighbours (not for the 'none' option).
seldist corresponds to the distance from the original control point when the 'none' option is not selected.

Inputs:

- filename: fullfile names of the files to be accessed in cell vector form. 
- ctlpoints: num_clt_points by 2 matrix with the coordinates of the 
     control points (pair of coordinated from where data will be extracted).
- DistClass: norm used to estimate the distance. 
	 This argument can take the next values:
		- Euclidean norm/distance: set {'Norm-2';'norm-2';2} and k=1 neighbour.
		- Infinite norm: set {inf;'inf'};
		- Minus infinite norm: set {-inf;'-inf'};
		- num-norm: for each integer num a norm can be defined. The case 2 
        is equivalent to the Euclidean norm or Euclidean Distance and k=1 neighbour.
- k: Number of neighbours (set k=[]; for 'none' option)
Example: * for 'none' option: [ncdata,selpoints]=ncextract(filename,ctlpoints,DistClass,k)
         * for other options: [ncdata,selpoints,distselind]=ncextract(filename,ctlpoints,DistClass,k)
Outputs:
- for 'none' option
     1. ncdata{j}(:,x) -> {j}: cell that corresponds to netcdf file currently
        accessed, (:,x): column of the time series extracted from the
        corresponding point
     2. selpoints -> num_ctl_points* by k matrix with the indexes of the actual control points.
- for all other options (k-nearest points)
    1. ncdata{j}(:,x,l) -> {j}: cell that corresponds to the netcdf file currently
        accessed, (:,x,l): column of the time series extracted from the
        corresponding k-nearest points to the point with selpoint{j}(x,:) row index x
    2. selpoints -> num_ctl_points*k by k matrix with the indexes of the k 
        nearest points to control points in ctl_points.
    3. distselind -> num_ctl_points by k matrix with the distance of the k
         nearest points to control points in ctl_points.

to do list:
 - automatic report generation
 - handle more dimensions like multiple grid levels of a netcdf file
 - sub-sample areas of the netcdf file
 - sub-sample time series of the netcdf file
 - mask NaN cells
 - improve processing speed via mexing or not ->pre-allocate all cells and matrixes
 - GUI

Coded by Voulanas Dimitrios on 2019-05-12
Copyright (c) 2019 Dimitrios Voulanas
http://https://github.com/dvoulanas
e-mail: dvoulanas@yahoo.com
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
  * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
  * The names of its contributor must not be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY Dimitrios Voulanas ''AS IS'' AND ANY EXPRESS OR 
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
IN NO EVENT SHALL Dimitrios Voulanas BE LIABLE FOR ANY DIRECT, INDIRECT, 
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT 
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
