# automated-nctools

automated-nctools is a suite of automated MATLAB functions. The first on is the one below.

1. function [ncdata,selpoints,distselind]=ncextract(filename,ctlpoints,DistClass,k)

This function searches for the nearest k-neighbors to the control points, within regular lonlat rectagular grid of any nc file or does
exact comparison to the grid cells of the gridded data
filename is a cell or string vector with fullfile names of the netcdf files to be accessed
ctlpoints are Num_ctl_points by 2 matrix. Each row it's a pair of lonlat coordinates.
layers_time is a 2 by 1 cell array if the nc file does not have a
layer/level dimension or a 2 by 2 cell array if the nc file has a
layer/level dimension
seldist is equal to the rows of ctlpoints
seldist corresponds to the degree distance from the original control point when the Distclass='none'; option
is not selected.
Also, this function features generation of process information that is printed at the command window.

You can add whatever variables you want, that the function produces, to the output of the function

Inputs:
 - filenames: fullfile names of the nc files to be accessed in cell or string vector form. 
 - ctlpoints: num_clt_points (total number of the control points pairs) by 2 matrix with the coordinates of the 
      control points (pair of coordinates from where data will be extracted, at [lon,lat] configuration).
 - DistClass: norm used to estimate the distance of k-nearest neighbors. 
	 This argument can take the next values:
		- Euclidean norm/distance: set {'Norm-2';'norm-2';2};
		- Infinite norm: set {inf;'inf'};
		- Minus infinite norm: set {-inf;'-inf'};
		- num-norm: for each integer num norm can be defined. The case 2 
         is equivalent to the Euclidean norm or Euclidean Distance.
       - k: Number of neighbors (set k=[]; for 'none' option)
Example: * for 'none' option: [ncdata,selpoints,sel_varNames,sel_dimNames,date_time,ncdimdata]=...
                              ncextract(filename,ctlpoints,DistClass,k)
         * for other options: [ncdata,selpoints,distselind,sel_varNames,sel_dimNames,date_time,ncdimdata]=...
                              ncextract(filename,ctlpoints,DistClass,k)
 Outputs:
 - for 'none' option
     1. ncdata{j}{x,:} -> {j}: cell that corresponds to netcdf file currently
        accessed, {x,:}: stored vector of the time series or matrix of the time series/layer extracted from the
         corresponding point (the times series is stored at each row while
         the layers/level is stored at the each column)
     2. selpoints -> num_ctl_points by ndim_ctl_points matrix with the indexes of the actual control points.
 - for all other options (k-nearest neigbour)
     1. ncdata{j}{x,l} -> {j}: cell that corresponds to the netcdf file currently
         accessed, {x,l}: stored vector of the time series or matrix of the time series/layer extracted from the
         corresponding k-nearest point to the point with selpoint{j}(x,:),
         row index x (the times series is stored at each row while
         the layers/level is stored at the each column)
     2. selpoints{j} -> num_ctl_points*k by num_ctl_points matrix with the indexes of the k-
         nearest points to control points in ctl_points, {j}: cell that corresponds to the netcdf file currently accessed.
	  3. distselind{j} -> num_ctl_points by k matrix with the distance of the
	  k-nearest points to control points in ctl_points, {j}: cell that corresponds to the netcdf file currently accessed.
     4. row{j} -> row indexes of the k-nearest neigbor at each netcdf, {j}: cell that corresponds to the netcdf file currently accessed.
     5. col{j} -> column indexes of the k-nearest neighor at each netcdf
     file, {j}: cell that corresponds to the netcdf file currently accessed.
 - for all options (can be skipped if unwanted): 
       1. sel_varNames{j} -> names of extracted variables -> {j}: cell that corresponds to the netcdf file currently
       2. sel_dimNames{j} -> names of extracted dimensions -> {j}: cell that corresponds to the netcdf file currently
       3. date_time{j} -> time vectors of each nc file in datetime format -> {j}: cell that corresponds to the netcdf file currently
       4. ncdimdata{j} -> dimension data of each nc file -> {j}: cell that corresponds to the netcdf file currently

 EXAMPLE (save to separate .m file)
 clc
 clear
 close all;
 ctlpoints=[12.9474,45.6161;12.36253738,45.0828757];
 % extract all layers/levels and time-series
 layers_time={1,1;Inf,Inf};
 % full paths of the file names
 filenames={'C:\yourpath...\filename.nc'};
 % euclidean norm
 DistClass=2;
 % 4 nearest neighbors
 k=4;
 % calls ncextract function
 [ncdata,selpoints,distselind,sel_varNames,sel_dimNames,date_time,ncdimdata]=...
     ncextract(filenames,ctlpoints,layers_time,DistClass,k);
 END EXAMPLE

 to do list:
 - automatic report generation
 - sub-sample areas of the netcdf file
 - add better date parser
 - mask NaN cells
 - extract user selected variable/s and not all variables
 - improve processing speed via mexing or not and pre-allocate all cells
   and matrixes
 - GUI

 v.2.0 Coded by Voulanas Dimitrios on 2020-03-25
 Copyright (c) 2020 Dimitrios Voulanas. 
 http://https://github.com/dvoulanas
 e-mail: dvoulanas@yahoo.com
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are 
 met:
 	* Redistributions of source code must retain the above copyright 
     notice, this list of conditions and the following disclaimer.
 	* Redistributions in binary form must reproduce the above copyright 
     notice, this list of conditions and the following disclaimer in the 
     documentation and/or other materials provided with the distribution.
 	* The names of its contributor must not be used to endorse or promote products 
     derived from this software without specific prior written permission.

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
