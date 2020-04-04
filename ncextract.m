function [ncdata,selpoints,sel_varNames,sel_dimNames,date_time,ncdimdata,distselind,row,col]=ncextract(filenames,ctlpoints,layers_time,DistClass,k)
%This function searches for the nearest k-neighbors to the control points, within regular lonlat rectagular grid of any nc file or does
%exact comparison to the grid cells of the gridded data
%filename is a cell or string vector with fullfile names of the netcdf files to be accessed
%ctlpoints are Num_ctl_points by 2 matrix. Each row it's a pair of lonlat coordinates.
%layers_time is a 2 by 1 cell array if the nc file does not have a
%layer/level dimension or a 2 by 2 cell array if the nc file has a
%layer/level dimension
%seldist is equal to the rows of ctlpoints
%seldist corresponds to the degree distance from the original control point when the Distclass='none'; option
%is not selected.
%Also, this function features generation of process information that is printed at the command
%window.
%
%You can add whatever variables you want, that the function produces, to the output
%of the function
%
%Inputs:
% - filenames: fullfile names of the nc files to be accessed in cell or string vector form. 
% - ctlpoints: num_clt_points (total number of the control points pairs) by 2 matrix with the coordinates of the 
%      control points (pair of coordinates from where data will be extracted, at [lon,lat] configuration).
% - DistClass: norm used to estimate the distance of k-nearest neighbors. 
%	 This argument can take the next values:
%		- Euclidean norm/distance: set {'Norm-2';'norm-2';2};
%		- Infinite norm: set {inf;'inf'};
%		- Minus infinite norm: set {-inf;'-inf'};
%		- num-norm: for each integer num norm can be defined. The case 2 
%         is equivalent to the Euclidean norm or Euclidean Distance.
%       - k: Number of neighbors (set k=[]; for 'none' option)
%Example: * for 'none' option: [ncdata,selpoints,sel_varNames,sel_dimNames,date_time,ncdimdata]=...
%                              ncextract(filename,ctlpoints,DistClass,k)
%         * for other options: [ncdata,selpoints,distselind,sel_varNames,sel_dimNames,date_time,ncdimdata]=...
%                              ncextract(filename,ctlpoints,DistClass,k)
% Outputs:
% - for 'none' option
%     1. ncdata{j}{x,:} -> {j}: cell that corresponds to netcdf file currently
%        accessed, {x,:}: stored vector of the time series or matrix of the time series/layer extracted from the
%         corresponding point (the times series is stored at each row while
%         the layers/level is stored at the each column)
%     2. selpoints -> num_ctl_points by ndim_ctl_points matrix with the indexes of the actual control points.
% - for all other options (k-nearest neigbour)
%     1. ncdata{j}{x,l} -> {j}: cell that corresponds to the netcdf file currently
%         accessed, {x,l}: stored vector of the time series or matrix of the time series/layer extracted from the
%         corresponding k-nearest point to the point with selpoint{j}(x,:),
%         row index x (the times series is stored at each row while
%         the layers/level is stored at the each column)
%     2. selpoints{j} -> num_ctl_points*k by num_ctl_points matrix with the indexes of the k-
%         nearest points to control points in ctl_points, {j}: cell that corresponds to the netcdf file currently accessed.
%	  3. distselind{j} -> num_ctl_points by k matrix with the distance of the
%	  k-nearest points to control points in ctl_points, {j}: cell that corresponds to the netcdf file currently accessed.
%     4. row{j} -> row indexes of the k-nearest neigbor at each netcdf, {j}: cell that corresponds to the netcdf file currently accessed.
%     5. col{j} -> column indexes of the k-nearest neighor at each netcdf
%     file, {j}: cell that corresponds to the netcdf file currently accessed.
% - for all options (can be skipped if unwanted): 
%       1. sel_varNames{j} -> names of extracted variables -> {j}: cell that corresponds to the netcdf file currently
%       2. sel_dimNames{j} -> names of extracted dimensions -> {j}: cell that corresponds to the netcdf file currently
%       3. date_time{j} -> time vectors of each nc file in datetime format -> {j}: cell that corresponds to the netcdf file currently
%       4. ncdimdata{j} -> dimension data of each nc file -> {j}: cell that corresponds to the netcdf file currently
%
% EXAMPLE (save to separate .m file)
% clc
% clear
% close all;
% ctlpoints=[12.9474,45.6161;12.36253738,45.0828757];
% % extract all layers/levels and time-series
% layers_time={1,1;Inf,Inf};
% % full paths of the file names
% filenames={'C:\yourpath...\filename.nc'};
% % euclidean norm
% DistClass=2;
% % 4 nearest neighbors
% k=4;
% % calls ncextract function
% [ncdata,selpoints,distselind,sel_varNames,sel_dimNames,date_time,ncdimdata]=...
%     ncextract(filenames,ctlpoints,layers_time,DistClass,k);
% END EXAMPLE
%
% to do list:
% - automatic report generation
% - sub-sample areas of the netcdf file
% - add better date parser
% - mask NaN cells
% - extract user selected variable/s and not all variables
% - improve processing speed via mexing or not and pre-allocate all cells
%   and matrixes
% - GUI
%
% v.2.0 Coded by Voulanas Dimitrios on 2020-03-25
% Copyright (c) 2020 Dimitrios Voulanas. 
% http://https://github.com/dvoulanas
% e-mail: dvoulanas@yahoo.com
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are 
% met:
% 	* Redistributions of source code must retain the above copyright 
%     notice, this list of conditions and the following disclaimer.
% 	* Redistributions in binary form must reproduce the above copyright 
%     notice, this list of conditions and the following disclaimer in the 
%     documentation and/or other materials provided with the distribution.
% 	* The names of its contributor must not be used to endorse or promote products 
%     derived from this software without specific prior written permission.

% THIS SOFTWARE IS PROVIDED BY Dimitrios Voulanas ''AS IS'' AND ANY EXPRESS OR 
% IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
% OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
% IN NO EVENT SHALL Dimitrios Voulanas BE LIABLE FOR ANY DIRECT, INDIRECT, 
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT 
% NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
% THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

timerValStart1=datetime(now,'ConvertFrom','datenum');
assert((iscell(filenames)|isstring(filenames))&isvector(filenames),'Filename must be cell or string vector.');
% %iterates through netcdf files
for j=1:size(filenames,2)
timerValStart=datetime(now,'ConvertFrom','datenum');
fprintf('%s: %i. Now processing %s.\n',datetime(now,'ConvertFrom','datenum'),j,filenames{j});  
assert(iscell(layers_time),'%s: %i. layers_time must be of cell type.\n',datetime(now,'ConvertFrom','datenum'),j);
finfo{j}=ncinfo(filenames{j});
ncvarNames{j}={finfo{j}.Variables.Name};
% %dim variables should only have 1 or 2 dimensions
% %data variables should have 3 or more dimensions
varNamelog{j}=cellfun(@numel,{finfo{j}.Variables.Size})>=3;
% %sets time as the first variable to be accessed
sel_dimNames{j}=sort(string(ncvarNames{j}(contains(ncvarNames{j},["lat","lon","depth","lvl","level","elevation"],'IgnoreCase',true))),'descend');
sel_varNames{j}=string(ncvarNames{j}(varNamelog{j}));
time_find{j}=contains(ncvarNames{j},["time"],'IgnoreCase',true) & cellfun(@numel,{finfo{j}.Variables.Size})==1;
sel_timeName{j}=string(ncvarNames{j}(time_find{j}));
fprintf('%s: %i. Now assessing dimension variables input.\n',datetime(now,'ConvertFrom','datenum'),j);
% % check if the nc file accessed has a layer/level dimension
if max(size(finfo{j}.Variables(varNamelog{j}).Size))-3>0
    assert(size(layers_time,1)==2 & size(layers_time,2)==2,...
        "This nc file has a layer/level dimension. layers_time vector must has rows: 2 for layer selection and cols: 2 for time selection.");
    assert((any(cellfun(@isdatetime,layers_time(:,2)))|any(cellfun(@isnumeric,layers_time(:,2))))&~any(cellfun(@isstring,layers_time(:,2))),...
        "Both elements of layers_time{:,2} must be datetime or numeric. Set layers_time{1,2} to a datetime between the start and end of the time-series of the gridded data, inclusive, or you can set it equal to 1 to begin extraction from the first value at the start of timeseries. Also, set layers_time{2,2} to a datetime between the value of layers_time{1,2} and end of the time-series of the gridded data, inclusive, or you can set it equal to Inf to extract data from the layers_time{1,2} value until the layers_time{2,2} value.");
    assert(all(cellfun(@isnumeric,layers_time(:,1))),...
        "Both elements of layers_time{:,1} must be numeric. Set layers_time{1,1} to a numeric between the first and last layer/level of the gridded data, inclusive, or you can set it equal to 1 to begin extraction from the first layer/level. Also, set layers_time{2,1} to a numeric value between the value of layers_time{1,2} and last layer/level of the gridded data, inclusive, or you can set it equal to Inf to extract data from the layers_time{1,1} layer/level until the last layer/level.");
    dim_check=1;
else
    assert(size(layers_time,1)==2 & size(layers_time,2)==1,...
        "This nc file has no layer/level dimension. layers_time vector must has rows: 2 and cols: 1 for time selection only.");
    assert(any(cellfun(@isdatetime,layers_time))|any(cellfun(@isnumeric,layers_time)),...
        "Both elements of layers_time must be datetime or numeric. Set layers_time{1,1} to a datetime between the start and end of the time-series of the gridded data, inclusive, or you can set it equal to 1 to begin extraction from the first value at the start of timeseries. Also, set layers_time{2,1} to a datetime between the value of layers_time{1,1} and end of the time-series of the gridded data, inclusive, or you can set it equal to Inf to extract data from the layers_time{2,1} value until the layers_time{2,1} value.");
    dim_check=0;
end
fprintf('%s: %i. Now reading time dimension: %s.\n',datetime(now,'ConvertFrom','datenum'),j,sel_timeName{j});
% % reads time dimension
temp=unique(ncread(filenames{j},sel_timeName{j}));
nctimedata{j}=sort(squeeze(temp(isfinite(temp)&~temp==0)),'ascend');
timeatt{j}={finfo{j}.Variables(time_find{j}).Attributes.Value};
timeatt_find=contains(timeatt{j},["since","from"],'IgnoreCase',true);
timeatt=split(string(timeatt{j}(find(timeatt_find,1)))," ");
timeatt_find=contains(timeatt,['-'],'IgnoreCase',true);
timeatt_date=strjoin(timeatt(find(timeatt_find):find(timeatt_find)+1));
%%%add what time format your file has here. There shouldn't be many cases
%%%around.
if string(timeatt(find(timeatt_find)+1))=="00:00:00"
    format='yyyy-MM-dd HH:mm:ss';
    fprintf('%s: %i. %s pattern has been found as the format of time variable.\n',datetime(now,'ConvertFrom','datenum'),j,format); 
elseif string(timeatt(find(timeatt_find)+1))=="00:00:00.0"
    format='yyyy-MM-dd HH:mm:ss.S';
    fprintf('%s: %i. %s pattern has been found as the format of time variable.\n',datetime(now,'ConvertFrom','datenum'),j,format); 
elseif string(timeatt(find(timeatt_find)+1))=="00:00"
    format='yyyy-MM-dd HH:mm';
    fprintf('%s: %i. %s pattern has been found as the format of time variable.\n',datetime(now,'ConvertFrom','datenum'),j,format); 
else
    format='yyyy-MM-dd HH:mm:ss.SSS';
    fprintf('%s: %i. %s pattern has been found as the format of time variable.\n',datetime(now,'ConvertFrom','datenum'),j,format); 
end
%%finds the time step of each dataset.
if string(timeatt{1})=="seconds"
    fprintf('%s: %i. %s has been found as variable time step.\n',datetime(now,'ConvertFrom','datenum'),j,"Seconds"); 
    date_time{j}=seconds(nctimedata{j})+datetime(timeatt_date,'Inputformat',format);
    date_time{j}.Format='dd-MMM-yyyy HH:mm:ss.SSS';
elseif string(timeatt{1})=="minutes"
    fprintf('%s: %i. %s has been found as variable time step.\n',datetime(now,'ConvertFrom','datenum'),j,"Minutes"); 
    date_time{j}=minutes(nctimedata{j})+datetime(timeatt_date,'Inputformat',format);
    date_time{j}.Format='dd-MMM-yyyy HH:mm:ss.SSS';
elseif string(timeatt{1})=="hours"
    fprintf('%s: %i. %s has been found as variable time step.\n',datetime(now,'ConvertFrom','datenum'),j,"Hours"); 
    date_time{j}=hours(nctimedata{j})+datetime(timeatt_date,'Inputformat',format);
    date_time{j}.Format='dd-MMM-yyyy HH:mm:ss.SSS';
elseif string(timeatt{1})=="days"
    fprintf('%s: %i. %s has been found as variable time step.\n',datetime(now,'ConvertFrom','datenum'),j,"Days"); 
    date_time{j}=caldays(nctimedata{j})+datetime(timeatt_date,'Inputformat',format);
    date_time{j}.Format='dd-MMM-yyyy HH:mm:ss.SSS';
elseif string(timeatt{1})=="months"
    fprintf('%s: %i. %s has been found as variable time step.\n',datetime(now,'ConvertFrom','datenum'),j,"Months"); 
    date_time{j}=calmonths(nctimedata{j})+datetime(timeatt_date,'Inputformat',format);
    date_time{j}.Format='MMM-yyyy';
else
    fprintf('%s: %i. %s has been found as variable time step.\n',datetime(now,'ConvertFrom','datenum'),j,"Years"); 
    date_time{j}=calyears(nctimedata{j})+datetime(timeatt_date,'Inputformat',format);
    date_time{j}.Format='yyyy';
end
% % converts datime to numeric indexes so the ncread function can access
% % the nc data for the time dimension
   layers_time(cellfun(@isdatetime,layers_time(:,size(layers_time,2))),size(layers_time,2))=...
       num2cell(find(contains(string(date_time{j}),string([layers_time{cellfun(@isdatetime,(layers_time(:,size(layers_time,2)))),size(layers_time,2)}]))));
   assert(all(cell2mat(layers_time(1,1:end))<=cell2mat(layers_time(2,1:end))),"Start extraction time and layer/level must less or equal to End extraction time and layer/level, respectively.");
% %reads dimension variables
for i=1:size(sel_dimNames{j},2)
   fprintf('%s: %i.%i. Now reading dimension: %s.\n',datetime(now,'ConvertFrom','datenum'),j,i,sel_dimNames{j}{i});
   temp1=unique(ncread(filenames{j},sel_dimNames{j}{i}));
   ncdimdata{j}{i}=sort(squeeze(temp1(isfinite(temp1)&~temp1==0)));
end
% %checks if case 'none'. This case searches for data in the corresponding individual grid
% %cells rather than the closest ones
if ~strcmp(DistClass,'none')
  [X,Y]=meshgrid(ncdimdata{j}{1},ncdimdata{j}{2});
  refpoints{j}=[reshape(X,1,[])',reshape(Y,1,[])'];
  [num_ctl_points,ndim_ctl_points]=size(ctlpoints);
  [num_ref_points,ndim_ref_points]=size(refpoints{j});
  distselind{j}=NaN(num_ctl_points,k);
  selind{j}=distselind{j};
  selpoints{j}=NaN(num_ctl_points*k,ndim_ctl_points);
  k_start=1;
  k_end=k;
else
 [num_ctl_points,ndim_ctl_points]=size(ctlpoints);
 selind{j}=NaN(num_ctl_points,1);
 selpoints{j}=NaN(num_ctl_points,ndim_ctl_points);
end
% %iterates through variables and control points
for f=1:size(sel_varNames{j},2)
    fprintf('%s: %i. Searching the indexes of selected grid points: %i in total.\n',datetime(now,'ConvertFrom','datenum'),j,size(ctlpoints,1)); 
    for x=1:size(selind{j},1)
        switch lower(DistClass)
          case {'none'}
            % %searches for the indexes of the grid points selected. If it doesnt find them it
            % %sets NaN as value        
            fprintf('%s: %i.%i. Matching coords XY: %10.5f %10.5f.\n',datetime(now,'ConvertFrom','datenum'),j,x,ctlpoints(x,:));
            pointLon=ncdimdata{j}{1}(ncdimdata{j}{1}>ctlpoints(x,1));
            if isempty(pointLon)
                selpoints{j}(x,1)=NaN;
                selind{j}(x,1)=NaN;
                fprintf('%s: %i.%i. Coord X:%10.5f is outside the grid.\n',datetime(now,'ConvertFrom','datenum'),j,x,ctlpoints(x,1));
            else
                selpoints{j}(x,1)=pointLon(1,1);
                selind{j}(x,1)=find(ncdimdata{j}{contains(sel_dimNames{j},["lon"],'IgnoreCase',true)}==selpoints{j}(x,1));
                fprintf('%s: %i.%i. Coord X:%10.5f has index %i.\n',datetime(now,'ConvertFrom','datenum'),j,x,selpoints{j}(x,1),selind{j}(x,1));
            end
            pointsLat=ncdimdata{j}{2}(ncdimdata{j}{2}>ctlpoints(x,2));
            if isempty(pointsLat)
                selpoints{j}(x,2)=NaN;
                selind{j}(x,2)=NaN;
                fprintf('%s: %i.%i. Coord Y: %10.5f is outside the grid.\n',datetime(now,'ConvertFrom','datenum'),j,x,ctlpoints(x,2));
            else
                selpoints{j}(x,2)=pointsLat(1,1);
                selind{j}(x,2)=find(ncdimdata{j}{2}==selpoints{j}(x,2));
                fprintf('%s: %i.%i. Coord Y: %10.5f has index %i.\n',datetime(now,'ConvertFrom','datenum'),j,x,selpoints{j}(x,2),selind{j}(x,2));
            end
            if isnan(selind{j}(x,1)) || isnan(selind{j}(x,2))
                fprintf('%s: %i.%i. Grid point XY: %10.5f %10.5f is outside the grid. Setting NaN.\n',datetime(now,'ConvertFrom','datenum'),j,x,ctlpoints(x,:));
                if dim_check==1
                    if layers_time{2,1}==Inf,layer_len=ncdimdata{j}{3}(end)-layers_time{1,1}+1;...
                    else,layer_len=layers_time{2,1}-layers_time{1,1}+1; end
                   ncdata{j}{f}{x,:}=NaN(size(nctimedata{j},1),layer_len);
                else
                   ncdata{j}{f}{x,:}=NaN(size(nctimedata{j},1),1);    
                end
            else
                fprintf('%s: %i.%i. Now extracting from variable %s at grid point XY:%10.5f%10.5f.\n',datetime(now,'ConvertFrom','datenum'),j,x,sel_varNames{j}{f},selpoints{j}(x,:));
                ncdata{j}{f}{x,:}=squeeze(ncread(filenames{j},sel_varNames{j}{f},[selind{j}(x,:),layers_time{1,:}],[1,1,layers_time{2,:}]))';      
            end
           case {'norm-2',2}
			 temp2=sqrt(sum((repmat(ctlpoints(x,:),num_ref_points,1)-refpoints{j}).^2,2));
		   case {inf,'inf'}
			 temp2=max(abs(repmat(ctlpoints(x,:),num_ref_points,1)-refpoints{j}),[],2);
		   case {-inf,'-inf'}
			 temp2=min(abs(repmat(ctlpoints(x,:),num_ref_points,1)-refpoints{j}),[],2);
           otherwise
	        if isnumeric(DistClass) && isfinite(DistClass)
			  temp2=sum(abs(repmat(ctlpoints(x,:),num_ref_points,1)-refpoints{j}).^DistClass,2).^(1/DistClass);
            else
		     error('Unknown norm: check the name or number of nearest neighbors.')
            end
        end
        if ~strcmp(DistClass,'none') && isfinite(k)
          [temp2,temp3]=sort(temp2);
          temp2=temp2(1:k)';
          temp3=temp3(1:k)';
          distselind{j}(x,:)=temp2;
          selind{j}(x,:)=find(~isnan(temp2));
        if ~isempty(selind{j}(x,:))
		  selind{j}(x,selind{j}(x,:))=temp3(selind{j}(x,:));
          selpoints{j}(k_start:k_end,:)=refpoints{j}(selind{j}(x,:),:);
          [row{j}(x,:),col{j}(x,:)]=ind2sub(size(X),selind{j}(x,:));
        end
        % % iterates through nearest neighbors
        for l=1:k
         if all(~isnan(row{j}(x,:)) & ~isnan(col{j}(x,:)))
           fprintf('%s: %i.%i.%i. Now extracting from variable %s at grid point XY:%10.5f%10.5f which has%10.5f deg distance from the original point XY:%10.5f%10.5f.\n',...
            datetime(now,'ConvertFrom','datenum'),j,x,l,sel_varNames{j}{f},selpoints{j}(k_start+l-1,:),distselind{j}(x,l),ctlpoints(x,:));     
           ncdata{j}{f}{x,l}=squeeze(ncread(filenames{j},sel_varNames{j}{f},[row{j}(x,l),col{j}(x,l),layers_time{1,:}],[1,1,layers_time{2,:}]))';    
         else
           fprintf('%s: %i.%i.%i. Grid point XY: %10.5f %10.5f is outside the grid. Setting NaN.\n',...
            datetime(now,'ConvertFrom','datenum'),j,x,l,selpoints{j}(k_start+l-1,:));
                if dim_check==1
                    if layers_time{2,1}==Inf,layer_len=ncdimdata{j}{3}(end)-layers_time{1,1}+1;...
                    else,layer_len=layers_time{2,1}-layers_time{1,1}+1; end
                   ncdata{j}{f}{x,l}=NaN(size(nctimedata{j},1),layer_len);
                else
                   ncdata{j}{f}{x,l}=NaN(size(nctimedata{j},1),1);    
                end
         end
        end
         k_start=k_start+k;
         k_end=k_end+k;
        end
    end
end
timerValEnd=datetime(now,'ConvertFrom','datenum');
dt=between(timerValStart,timerValEnd);
fprintf('%s: %i. Extraction for file: %s was completed on %s, after %s of elapsed time.\n',datetime(now,'ConvertFrom','datenum'),j,filenames{j},timerValEnd,dt);
end
timerValEnd1=datetime(now,'ConvertFrom','datenum');
dt1=between(timerValStart1,timerValEnd1);
fprintf('%s: Extraction was completed after %s in total.\n',datetime(now,'ConvertFrom','datenum'),dt1);
end
