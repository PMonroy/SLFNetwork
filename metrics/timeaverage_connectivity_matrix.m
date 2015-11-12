clear all
close all

datestr(now)

addpath /home/pmonroy/MATLAB/function/


output_dir='/home/pmonroy/ESCOLA/RUNS/SLFNetwork/';

% Load bathymetry
bathymetry=load('bathy_BENGUELA_etopo1.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Users inputs for which effect to test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%networkdomain1_vflow1_depth1000_nodesize0167.grid
%networkdomain1_vflow1_startdepth0025_finaldepth1000_particlespacing00100_nodesize0167_startdate01030008_vsink050_intstep0025.matrix

% Manual inputs: parameters
% Low sinking
%d='1';r='0167';p='00100';it='0025';zinit='0025';zend='1000';vsink='010';
% Moderate sinking
%d='1';r='0167';p='00100';it='0025';zinit='0025';zend='1000';vsink='050';
% Fast sinking
d='1';r='0167';p='00100';it='0025';zinit='0025';zend='1000';vsink='100';

% Upwelling relaxation season, both years together
%year={'0008' '0009'};daymonth={'0103' '1003' '2003' '3003' '1004' '2004' '3004'};
%figname = ['SLFN_vsink',vsink,'_Mean_relaxation_both_years.pdf'];

% Upwelling season, both years together
%year={'0008' '0009'};daymonth={'0109' '1009' '2009' '3009' '1010' '2010' '3010'};
%figname = ['SLFN_vsink',vsink,'_Mean_upwelling_both_years.pdf'];

% Upwelling relaxation season, only 1 year
%year={'0008'};daymonth={'0103' '1003' '2003' '3003' '1004' '2004' '3004'};
%figname = ['SLFN_vsink',vsink,'_Mean_relaxation_one_year.pdf'];

% Upwelling season, only 1 year
year={'0008'};daymonth={'0109' '1009' '2009' '3009' '1010' '2010' '3010'};
figname = ['SLFN_vsink',vsink,'_Mean_upwelling_one_year.pdf'];

% Build dt from year and t0
for i=1:length(year)
    if i==1
        startdate=strcat(daymonth,year(i));
    else
        startdate=cat(2,startdate,strcat(daymonth,year(i)));
    end;
end;


				% CONSTRUCT file names and LOAD data

				% Load grid nodes
sfilegrid=['networkdomain',d,'_vflow1_depth',zend,'_nodesize',r,'.grid'];
grid=load([output_dir,sfilegrid]);


				% CONSTRUCT file names and LOAD data
 cmatrix.npoint=zeros([length(grid),length(grid),length(startdate)]);
 cmatrix.meantime=zeros([length(grid),length(grid),length(startdate)]);
 cmatrix.vartime=zeros([length(grid),length(grid),length(startdate)]);
 cmatrix.node=zeros([length(grid),length(grid),length(startdate)]);

 for i=1:length(startdate)
				% Building files names
   sfilematrix{i}=['networkdomain',d,'_vflow1_startdepth',zinit,'_finaldepth',zend,'_particlespacing',p,'_nodesize',r,'_startdate',startdate{i},'_vsink',vsink,'_intstep',it,'.matrix'];
				% Load data
   listlinks=load([output_dir,sfilematrix{i}]);
				% Make Connectivity list Matlab friendly
   listlinks(:,1:2)=listlinks(:,1:2)+1; 
   
   for j=1:length(listlinks)
     startnode=listlinks(j,1);
     finalnode=listlinks(j,2);

     cmatrix.npoint(startnode,finalnode,i) = listlinks(j,3);
     cmatrix.meantime(startnode,finalnode,i) = listlinks(j,4);
     cmatrix.vartime(startnode,finalnode,i) = listlinks(j,5);
     cmatrix.node(startnode,finalnode,i) = 1;
   end
   dispprocess(i,length(startdate));
   clear listlinks;
 end

 Meancmatrix.npoint=sum(cmatrix.npoint,3)/length(startdate);
 Meancmatrix.node=sum(cmatrix.node,3)/length(startdate);

 Instrength.npoint=sum(Meancmatrix.npoint,1);
 Outstrength.npoint=sum(Meancmatrix.npoint,2);

 Indegree=sum(Meancmatrix.node,1);
 Outdegree=sum(Meancmatrix.node,2);

 % MEAN TIME = PUFFFF
 %Meancmatrix.meantime=sum(cmatrix.meantime,3)/length(startdate);
 %Meancmatrix.vartime=sum(cmatrix.vartime,3)/length(startdate);
 %Instrength.meantime=sum(Meancmatrix.meantime,1);
 %Outstrength.meantime=sum(Meancmatrix.meantime,2);
 %Instrength.vartime=sum(Meancmatrix.vartime,1);
 %Outstrength.vartime=sum(Meancmatrix.vartime,2);

 datestr(now)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;

subplot(2,2,1);
hold on;
deltalat = (str2double(r)/1000);
rectangles = [grid(:,1)-(grid(:,3)/2), grid(:,2)-(deltalat/2), grid(:,3), deltalat*ones(size(grid,1),1)];

rectanglecolor =  colormap(jet);
m = length(rectanglecolor);
cmin = min(Indegree(1,:));
cmax = max(Indegree(1,:));
for i=1:length(rectangles)
  index = fix((Indegree(1,i)-cmin)/(cmax-cmin)*m)+1;
				%Clamp values outside the range [1 m]
  index(index<1) = 1;
  index(index>m) = m; 
  rectangle('Position', rectangles(i,:),'FaceColor', rectanglecolor(index,:));
end;
caxis([cmin cmax]);
colormap(jet);
%adds the colorbar to your plot
colorbar

title(['Indegree'],'fontsize',10);
xlabel('Longitude\circE','fontsize',8);
ylabel('Latitude\circN','fontsize',8);
set(gca,'xlim',[min(grid(:,1)) max(grid(:,1))],'ylim',[min(grid(:,2)) max(grid(:,2))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,2);
hold on;
deltalat = (str2double(r)/1000);
rectangles = [grid(:,1)-(grid(:,3)/2), grid(:,2)-(deltalat/2), grid(:,3), deltalat*ones(size(grid,1),1)];

rectanglecolor =  colormap(jet);
m = length(rectanglecolor);
cmin = min(Outdegree(:,1));
cmax = max(Outdegree(:,1));
for i=1:length(rectangles)
    index = fix((Outdegree(i,1)-cmin)/(cmax-cmin)*m)+1;
    %Clamp values outside the range [1 m]
    index(index<1) = 1;
    index(index>m) = m; 
    
    rectangle('Position', rectangles(i,:),'FaceColor', rectanglecolor(index,:));
end;
caxis([cmin cmax]);
colormap(jet);
%adds the colorbar to your plot
colorbar

title(['Outdegree'],'fontsize',10);
xlabel('Longitude\circE','fontsize',8);
ylabel('Latitude\circN','fontsize',8);
set(gca,'xlim',[min(grid(:,1)) max(grid(:,1))],'ylim',[min(grid(:,2)) max(grid(:,2))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3);
hold on;
deltalat = (str2double(r)/1000);
rectangles = [grid(:,1)-(grid(:,3)/2), grid(:,2)-(deltalat/2), grid(:,3), deltalat*ones(size(grid,1),1)];

rectanglecolor =  colormap(jet);
m = length(rectanglecolor);
cmin = min(Instrength.npoint(1,:));
cmax = max(Instrength.npoint(1,:));
for i=1:length(rectangles)
  index = fix((Instrength.npoint(1,i)-cmin)/(cmax-cmin)*m)+1;
				%Clamp values outside the range [1 m]
  index(index<1) = 1;
  index(index>m) = m; 
  rectangle('Position', rectangles(i,:),'FaceColor', rectanglecolor(index,:));
end;
caxis([cmin cmax]);
colormap(jet);
%adds the colorbar to your plot
colorbar

title(['Instrength'],'fontsize',10);
xlabel('Longitude\circE','fontsize',8);
ylabel('Latitude\circN','fontsize',8);
set(gca,'xlim',[min(grid(:,1)) max(grid(:,1))],'ylim',[min(grid(:,2)) max(grid(:,2))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,4);
hold on;
deltalat = (str2double(r)/1000);
rectangles = [grid(:,1)-(grid(:,3)/2), grid(:,2)-(deltalat/2), grid(:,3), deltalat*ones(size(grid,1),1)];

rectanglecolor =  colormap(jet);
m = length(rectanglecolor);
cmin = min(Outstrength.npoint(:,1));
cmax = max(Outstrength.npoint(:,1));
for i=1:length(rectangles)
    index = fix((Outstrength.npoint(i,1)-cmin)/(cmax-cmin)*m)+1;
    %Clamp values outside the range [1 m]
    index(index<1) = 1;
    index(index>m) = m; 
    
    rectangle('Position', rectangles(i,:),'FaceColor', rectanglecolor(index,:));
end;
caxis([cmin cmax]);
colormap(jet);
%adds the colorbar to your plot
colorbar

title(['Outstrength'],'fontsize',10);
xlabel('Longitude\circE','fontsize',8);
ylabel('Latitude\circN','fontsize',8);
set(gca,'xlim',[min(grid(:,1)) max(grid(:,1))],'ylim',[min(grid(:,2)) max(grid(:,2))]);

set(gcf,'PaperOrientation','landscape')
set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperType','A4')
set(gcf,'PaperPosition',[.5 .5 28. 18.])

print('-dpdf',figname);