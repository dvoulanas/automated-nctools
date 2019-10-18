function [ncdata,selpoints,distselind]=ncextract(filename,ctlpoints,DistClass,k)
%This function searches for the nearest k-neighbors to the control points, within XY or does
%exact comparison to the grid cells of the gridded data
%coords of regular lonlat rectagular grid of netcdfs files or makes exact comparisons
%for each corresponding grid cell
%filename is a cell vector with fullfile name of the netcdf fiel to be accessed
%ctlpoints are Num_ctl_points by 2 matrix. Each row it's a pair of XY coordinates.
%seldist is equal to the rows of ctlpoints, and the number of columns is
%k, which corresponds to the number of neighbours (not for the 'none' option).
%seldist corresponds to the distance from the original control point when the 'none' option
%is not selected.
%
%Inputs:
% - filename: fullfile names of the files to be accessed in cell vector form. 
% - ctlpoints: num_clt_points by 2 matrix with the coordinates of the 
%      control points (pair of coordinated from where data will be extracted).
% - DistClass: norm used to estimate the distance. 
%	 This argument can take the next values:
%		- Euclidean norm/distance: set {'Norm-2';'norm-2';2} and k=1 neighbour.
%		- Infinite norm: set {inf;'inf'};
%		- Minus infinite norm: set {-inf;'-inf'};
%		- num-norm: for each integer num a norm can be defined. The case 2 
%         is equivalent to the Euclidean norm or Euclidean Distance and k=1 neighbour.
% - k: Number of neighbours (set k=[]; for 'none' option)
%Example: * for 'none' option: [ncdata,selpoints]=ncextract(filename,ctlpoints,DistClass,k)
%         * for other options: [ncdata,selpoints,distselind]=ncextract(filename,ctlpoints,DistClass,k)
% Outputs:
% - for 'none' option
%     1. ncdata{j}(:,x) -> {j}: cell that corresponds to netcdf file currently
%        accessed, (:,x): column of the time series extracted from the
%        corresponding point
%     2. selpoints -> num_ctl_points* by k matrix with the indexes of the actual control points.
% - for all other options (k-nearest points)
%     1. ncdata{j}(:,x,l) -> {j}: cell that corresponds to the netcdf file currently
%         accessed, (:,x,l): column of the time series extracted from the
%         corresponding k-nearest points to the point with selpoint{j}(x,:) row index x
%     2. selpoints -> num_ctl_points*k by k matrix with the indexes of the k 
%         nearest points to control points in ctl_points.
%	  3. distselind -> num_ctl_points by k matrix with the distance of the k
%         nearest points to control points in ctl_points.

% to do list:
% - automatic report generation
% - handle more dimensions like multiple grid levels of a netcdf file
% - sub-sample areas of the netcdf file
% - sub-sample time series of the netcdf file
% - mask NaN cells
% - improve processing speed via mexing or not ->pre-allocate all cells
% and matrixes
% - GUI
%
% Coded by Voulanas Dimitrios on 2019-05-12
% Copyright (c) 2019 Dimitrios Voulanas. 
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

% %iterates through netcdf files
for j=1:size(filename,2)
fprintf('%s: %i. Now processing %s.\n',datetime(now,'ConvertFrom','datenum'),j,filename{j});  
finfo{j}=ncinfo(filename{j});
ncvarNames{j}={finfo{j}.Variables.Name};
% %dim variable should only have 1 dimension
ncdimlog{j}=cellfun(@numel,{finfo{j}.Variables.Size})==1;
% %data variables should have 3 or more dimensions
varNamelog{j}=cellfun(@numel,{finfo{j}.Variables.Size})>=3;
% %sets time as the first variable to be accessed
sel_dimNames{j}=sort(string(ncvarNames{j}(ncdimlog{j})),'descend');
sel_varNames{j}=string(ncvarNames{j}(varNamelog{j}));
% %reads dimension variable as set above
for i=1:size(sel_dimNames{j},2)
   fprintf('%s: %i. Now reading dimension variable: %s.\n',datetime(now,'ConvertFrom','datenum'),i,sel_dimNames{j}{i}); 
   ncdimdata{j}{i}=sort(squeeze(ncread(filename{j},sel_dimNames{j}{i})));
end
timeatt_full{j}={finfo{j}.Variables(1).Attributes.Value};
timeatt_find_full=strfind(timeatt_full{j},'-');
emptyCells_timeatt_full=~cellfun(@isempty,timeatt_find_full);
timeatt=split(char(timeatt_full{j}(emptyCells_timeatt_full)));
timeatt_find=strfind(timeatt,'-');
emptyCells_timeatt=~cellfun(@isempty,timeatt_find);
timeatt_date=strjoin(timeatt(find(emptyCells_timeatt):find(emptyCells_timeatt)+1));
%%%add what time format your file has here. There shouldn't be many cases
%%%around.
if string(timeatt(find(emptyCells_timeatt)+1))=="00:00:00"
    format='yyyy-MM-dd HH:mm:ss';
    fprintf('%s: %s pattern has been found as the format of time variable.\n',datetime(now,'ConvertFrom','datenum'), format); 
elseif string(timeatt(find(emptyCells_timeatt)+1))=="00:00:00.0"
    format='yyyy-MM-dd HH:mm:ss.S';
    fprintf('%s: %s pattern has been found as the format of time variable.\n',datetime(now,'ConvertFrom','datenum'),format); 
elseif string(timeatt(find(emptyCells_timeatt)+1))=="00:00"
    format='yyyy-MM-dd HH:mm';
    fprintf('%s: %s pattern has been found as the format of time variable.\n',datetime(now,'ConvertFrom','datenum'), format); 
else
    format='yyyy-MM-dd HH:mm:ss.SSS';
    fprintf('%s: %s pattern has been found as the format of time variable.\n',datetime(now,'ConvertFrom','datenum'),format); 
end
%%finds the time step of dataset.
if string(timeatt{1})=="seconds"
    fprintf('%s: %s has been found as variable time step.\n',datetime(now,'ConvertFrom','datenum'),"Seconds"); 
    time_datenum{j}=hours(ncdimdata{j}{1})+datetime(timeatt_date,'Inputformat',format);
    date_time_datenum{j}{i}.Format = 'dd-MMM-yyyy HH:mm:ss.SSS';
elseif string(timeatt{1})=="hours"
    fprintf('%s: %s has been found as variable time step.\n',datetime(now,'ConvertFrom','datenum'), "Hours"); 
    time_datenum{j}=hours(ncdimdata{j}{1})+datetime(timeatt_date,'Inputformat',format);
    date_time_datenum{j}{i}.Format = 'dd-MMM-yyyy HH:mm:ss.SSS';
elseif string(timeatt{1})=="days"
    fprintf('%s: %s has been found as variable time step.\n',datetime(now,'ConvertFrom','datenum'),"Days"); 
    time_datenum{j}=caldays(ncdimdata{j}{1})+datetime(timeatt_date,'Inputformat',format);
    date_time_datenum{j}{i}.Format = 'dd-MMM-yyyy HH:mm:ss.SSS';
elseif string(timeatt{1})=="months"
    fprintf('%s: %s has been found as variable time step.\n',datetime(now,'ConvertFrom','datenum'), "Months"); 
    time_datenum{j}=calmonths(ncdimdata{j}{1})+datetime(timeatt_date,'Inputformat',format);
    date_time_datenum{j}{i}.Format = 'MMM-yyyy';
else
    fprintf('%s: %s has been found as variable time step.\n',datetime(now,'ConvertFrom','datenum'), "Years"); 
    time_datenum{j}=calyears(ncdimdata{j}{1})+datetime(timeatt_date,'Inputformat',format);
    date_time_datenum{j}{i}.Format = 'yyyy';
end
% %check if case 'none'. This case searches for data in the corresponding individual grid
% %cells rather than the closest ones
if ~strcmp(DistClass,'none')
  [X,Y]=meshgrid(ncdimdata{1}{2},ncdimdata{1}{3});
  refpoints{j}=[reshape(X,1,[])',reshape(Y,1,[])'];
  [num_ctl_points,ndim_ctl_points]=size(ctlpoints);
  [num_ref_points,ndim_ref_points]=size(refpoints{j});
  distselind{j}=NaN(num_ctl_points,k);
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
    fprintf('%s: Searching the indexes of selected grid points: %i in total.\n',datetime(now,'ConvertFrom','datenum'),size(ctlpoints,1)); 
    for x=1:size(selind{j},1)
        switch lower(DistClass)
          case {'none'}
            % %searches for the indexes of the grid points selected. If it doesnt find them it
            % %sets NaN as value        
            fprintf('%s: %i. Matching coords XY: %10.5f %10.5f.\n',datetime(now,'ConvertFrom','datenum'),x,ctlpoints(x,:));
            pointLon=ncdimdata{j}{2}(ncdimdata{j}{2}>ctlpoints(x,1));
            if isempty(pointLon)
                selpoints{j}(x,1)=NaN;
                selind{j}(x,1)=NaN;
                fprintf('%s: %i. Coord X:%10.5f is outside the grid.\n',datetime(now,'ConvertFrom','datenum'),x,ctlpoints(x,1));
            else
                selpoints{j}(x,1)=pointLon(1,1);
                selind{j}(x,1)=find(ncdimdata{j}{2}==selpoints{j}(x,1));
                fprintf('%s: %i. Coord X:%10.5f has index %i.\n',datetime(now,'ConvertFrom','datenum'),x,selpoints{j}(x,1),selind{j}(x,1));
            end
            pointsLat=ncdimdata{j}{3}(ncdimdata{j}{3}>ctlpoints(x,2));
            if isempty(pointsLat)
                selpoints{j}(x,2)=NaN;
                selind{j}(x,2)=NaN;
                fprintf('%s: %i. Coord Y: %10.5f is outside the grid.\n',datetime(now,'ConvertFrom','datenum'),x,ctlpoints(x,2));
            else
                selpoints{j}(x,2)=pointsLat(1,1);
                selind{j}(x,2)=find(ncdimdata{j}{3}==selpoints{j}(x,2));
                fprintf('%s: %i. Coord Y: %10.5f has index %i.\n',datetime(now,'ConvertFrom','datenum'),x,selpoints{j}(x,2),selind{j}(x,2));
            end
            if isnan(selind{j}(x,1)) || isnan(selind{j}(x,2))
                fprintf('%s: %i. Grid point XY: %10.5f %10.5f is outside the grid. Setting NaN.\n',datetime(now,'ConvertFrom','datenum'),x,ctlpoints(x,:));
                ncdata{j}{f}(:,x)=NaN(size(ncdimdata{j}{1},1),1); 
            else
                fprintf('%s: %i. Now extracting from variable %s at grid point XY:%10.5f%10.5f.\n',datetime(now,'ConvertFrom','datenum'),x,sel_varNames{j}{f},selpoints{j}(x,:));
                ncdata{j}{f}(:,x)=squeeze(ncread(filename{j},sel_varNames{j}{f},[selind{j}(x,:),1],[1,1,Inf]));      
            end
           case {'norm-2',2}
			 temp=sqrt(sum((repmat(ctlpoints(x,:),num_ref_points,1)-refpoints{j}).^2,2));
		   case {inf,'inf'}
			 temp=max(abs(repmat(ctlpoints(x,:),Nest_ref_points,1)-refpoints{j}),[],2);
		   case {-inf,'-inf'}
			 temp=min(abs(repmat(ctlpoints(x,:),Nest_ref_points,1)-refpoints{j}),[],2);
           otherwise
	       if isnumeric(DistClass)
			  temp=sum(abs(repmat(ctlpoints(x,:),Nest_ref_points,1)-refpoints{j}).^DistClass,2).^(1/DistClass);
           else
		     error('Unknown norm: check the name')
           end
        end
        if ~strcmp(DistClass,'none')
          [temp,temp1]=sort(temp);
          temp=temp(1:k)';
          temp1=temp1(1:k)';
          distselind{j}(x,:)=temp;
          selind{j}(x,:)=find(~isnan(temp));
        if ~isempty(selind{j}(x,:))
		  selind{j}(x,selind{j}(x,:))=temp1(selind{j}(x,:));
          selpoints{j}(k_start:k_end,:)=refpoints{j}(selind{j}(x,:),:);
          [row{j}(x,:),col{j}(x,:)]=ind2sub(size(X),selind{j}(x,:));
        end
        for l=1:k
         if all(~isnan(row{j}(x,:)) & ~isnan(col{j}(x,:)))
           fprintf('%s: %i.%i. Now extracting from variable %s at grid point XY:%10.5f%10.5f which has%10.5f deg distance from the original point XY:%10.5f%10.5f.\n',...
            datetime(now,'ConvertFrom','datenum'),x,l,sel_varNames{j}{f},selpoints{j}(k_start+l-1,:),distselind{j}(x,l),ctlpoints(x,:));     
           ncdata{j}{f}(:,x,l)=squeeze(ncread(filename{j},sel_varNames{j}{f},[row{j}(x,l),col{j}(x,l),1],[1,1,Inf]));    
         else
           fprintf('%s: %i.%i. Grid point XY: %10.5f %10.5f is outside the grid. Setting NaN.\n',...
            datetime(now,'ConvertFrom','datenum'),x,l,selpoints{j}(k_start+l-1,:));
           ncdata{j}{f}(:,x)=NaN;
         end
        end
         k_start=k_start+k;
         k_end=k_end+k;
        end
    end
end
end
end