% created by JN Kather and CA Weis 2015-2016
% jakob.kather@medma.uni-heidelberg.de
% license: MIT license, see separate file for license and disclaimer
%
% parts of the code are based on the following publication
% Kather, JN et al. Continuous representation of tumor microvessel 
% density and detection of angiogenic hotspots in histological 
% whole-slide images. Oncotarget 5, (2015). DOI: 10.18632/oncotarget.4383
%
% others are based on this publication:
% Kather, JN et al. New Colors for Histology: Optimized Bivariate 
% Color Maps Increase Perceptual Contrast in Histological Images. 
% PLoS One 10, e0145572 (2015). DOI: 10.1371/journal.pone.0145572
%
% File description:
% this program loads hotspot data and vessel coordinates and computes
% statistics such as mean MVD per hotspot and hotspots per area. It has to
% be run manually.
%%
% -- initialize program
clear all, close all; 
format compact; clc; 

% -- add toolboxes to current path
addpath([pwd,'/subroutines'],'-end'); % my own subroutines
addpath([pwd,'/poly_stuff'],'-end');  % code from http://www.mathworks.com/matlabcentral/fileexchange/10391-fast-points-in-polygon-test
addpath([pwd,'/Yen'],'-end');  % code from http://www.mathworks.com/matlabcentral/fileexchange/10391-fast-points-in-polygon-test

cnst = getConstants(); % define constants
cnst.path = strcat(cnst.img_folder);
printFigs = false; % save figures?
showResult = false; % show figures?

% specify directories and file to load, read all files
%id = '*';
filtr = '*'; % default *
allFiles = dir([cnst.cachepath,'*',filtr,'*_vessel_measurements.mat']);
first = 1;

% do this for every suitable file
for i_Files=1:numel(allFiles)
    
currMatFile = allFiles(i_Files).name; 

disp(['#######################',10,currMatFile,10]);
load([cnst.cachepath,currMatFile]); % load general data

try
    load([cnst.cachepath,currImage.tiffname,'_ROI_data.mat']); % load ROI data
catch
    warning('could not load ROI data. break');
    continue
end

warning('loaded saved variables from disk');

% iterate through all ROIs in image
for i_ROI=1:numel(ROINameContainer)
    disp(['ROI --> ', char(ROINameContainer{i_ROI})]);
    
    if strcmp(char(ROINameContainer{i_ROI}),'FullTumorValidated')
        
    currPoly = polygonContainer{i_ROI};
    
    % load hotspot map data
    disp('loading hotspot data for current ROI');
    load([cnst.cachepath,currImage.tiffname,'_hotspot_measurements_',...
        char(ROINameContainer{i_ROI}),'.mat']);  % load hotspot data 
    
    % extract in-polygon coordinates of blood vessels
    VessInP = inpolygon(currImage.ObjTable.Centroid(:,2)/cnst.resizefactor,...
                    currImage.ObjTable.Centroid(:,1)/cnst.resizefactor,...
                    currPoly(:,2),currPoly(:,1));
    
    if showResult
    h1 = figure();
    subplot(1,2,1)
    patch(currPoly(:,1),...
          currPoly(:,2),'red',...
          'FaceAlpha',0.1,'EdgeColor','black','LineWidth',2);    
    hold on
     
%     contourf(currMap.X,currMap.Y,currMap.densityOBS,...
%         cnst.contourlines,'LineStyle','none');
%     colormap parula
    hold on
    
    scatter(currImage.ObjTable.Centroid(VessInP,1)/cnst.resizefactor,...
        currImage.ObjTable.Centroid(VessInP,2)/cnst.resizefactor,2,'black','filled');
    title('objects in ROI');
%     scatter(currImage.ObjTable.Centroid(:,1)/cnst.resizefactor,...
%         currImage.ObjTable.Centroid(:,2)/cnst.resizefactor,20,'blue','filled');
    
    axis equal tight on;
    view(0,90);  
    
    subplot(1,2,2)
    contourf(currMap.X,currMap.Y,currMap.densityOBS,...
        cnst.contourlines,'LineStyle','none');
    title('vessel density');
    colormap hot
    caxis([0 4e-08])
    colorbar
      axis equal tight on;
    view(0,90);  
    
    set(gcf,'Color','w');
    
    suptitle(strrep([currImage.tiffname,10,char(ROINameContainer{i_ROI})],'_',' '));
    drawnow
    
    if printFigs
    print(h1,'-dpng','-r450',[cnst.resultsFolder,currImage.tiffname,...
           char(ROINameContainer{i_ROI}), '_05_OVERVIEW_', '.png']);  disp('done printing');
    end
    end
    
    % save quantitative results to table
    [currentTable, h13] = getHotspotStatistics(currImage.ObjTable.Centroid(:,1)/cnst.resizefactor,...
                                             currImage.ObjTable.Centroid(:,2)/cnst.resizefactor,...
                                             currPoly(:,1),...
                                             currPoly(:,2),...
                                             currMap.X,...
                                             currMap.Y,...
                                             currMap.densityOBS,...
                                             currMap.probabilityMap,...
                                             currMap.sigLevelBonferroni,...
                                             currImage.tiffname,...
                                             char(ROINameContainer{i_ROI}),...
                                             currImage.ObjTable, cnst);
    
    % --- get cumulative positive excess of distance histograms - VESSELS
    histogrTitle = 'Lumen';
    try
        load([cnst.cachepath,currImage.tiffname,...
                '_distances-histogram-class_',histogrTitle,'.mat']); 
    catch
        error(['did not find histogram file: ',cnst.cachepath,currImage.tiffname,...
                '_distances-histogram-class_',histogrTitle,'.mat']); 
    end
    
    try
        currentTable.LumenZoneWidth = meanZoneWidth;
        currentTable.LumenCumulExcess = AverageCumulExcess;
        currentTable.LumenExcessLoCI = cumulExcess_loCI;
        currentTable.LumenExcessHiCI = cumulExcess_hiCI;
        currentTable.LumenExcessStd = CumulExcessStd;
		currentTable.LumenProbValue = ttestPVal;
        currentTable.LumenProbDecis = ttestDecis;
        currentTable.LumenExcessDirection = meanExcessDirection;
        currentTable.LumenExists = exist_Lumen;
    catch
        warning('could not assign values, set 0 or -1');
        currentTable.LumenZoneWidth = -1;
        currentTable.LumenCumulExcess = 0;
		currentTable.LumenProbValue = -1;
        currentTable.LumenExcessLoCI = 0;
        currentTable.LumenExcessHiCI = 0;
        currentTable.LumenExcessStd = 0;
		currentTable.LumenProbValue = -1;
        currentTable.LumenProbDecis = -1;
        currentTable.LumenExcessDirection = 0;
        currentTable.LumenExists = -1;     
        pause(3);
    end
    
    clear ExcessDirection cumulExcess zoneWidth;
    
    % --
    histogrTitle = 'Adjacent Tissue';
    try
    load([cnst.cachepath,currImage.tiffname,...
            '_distances-histogram-class_',histogrTitle,'.mat']); 
    catch
        error('did not find histogram file');
    end    
    
    try
        currentTable.AdjacentZoneWidth = meanZoneWidth;
        currentTable.AdjacentCumulExcess = AverageCumulExcess;
        currentTable.AdjacentExcessLoCI = cumulExcess_loCI;
        currentTable.AdjacentExcessHiCI = cumulExcess_hiCI;
        currentTable.AdjacentExcessStd = CumulExcessStd;
		currentTable.AdjacentProbValue = ttestPVal;
        currentTable.AdjacentProbDecis = ttestDecis;
        currentTable.AdjacentExcessDirection = meanExcessDirection;
        currentTable.AdjacentExists = exist_Adjacent;
    catch
        warning('could not assign values, set 0 or -1');
        currentTable.AdjacentZoneWidth = -1;
        currentTable.AdjacentCumulExcess = 0;
		currentTable.AdjacentProbValue = -1;
        currentTable.AdjacentExcessLoCI = 0;
        currentTable.AdjacentExcessHiCI = 0;
        currentTable.AdjacentExcessStd = 0;
		currentTable.AdjacentProbValue = -1;
        currentTable.AdjacentProbDecis = -1;
        currentTable.AdjacentExcessDirection = 0;
        currentTable.AdjacentExists = -1;     
        pause(3);
    end
        
    clear ExcessDirection cumulExcess exist_Lumen exist_Adjacent zoneWidth;
    % ---         
   
    if printFigs
    print(h13,'-dpng','-r450',[cnst.resultsFolder,currImage.tiffname,...
           char(ROINameContainer{i_ROI}), '_06_HOTSPOTS_', '.png']);  disp('done printing');
    end
     
    %----------------
    % load CD03 and Ki67 data
    disp('loading CD3/Ki67 data...');
    
    for currDataName = {'CD03_DATA','Ki67_DATA'}
        for histogrTitle = {'Lumen','Adjacent_Tissue'}
            try
                disp(['trying to load ... ',10, cnst.cachepath,currImage.tiffname,...
                '_distances-histogram-Data',char(currDataName),...
                'class_',strrep(char(histogrTitle),'_',' '),'.mat']);
            
                load([cnst.cachepath,currImage.tiffname,...
                '_distances-histogram-Data',char(currDataName),...
                'class_',strrep(char(histogrTitle),'_',' '),'.mat']);
        
                % save to table
                
                %cumulExcess, zoneWidth, ExcessDirection, exist_Lumen, exist_Adjacent
                % clean up
                currentTable.([char(currDataName),'_',char(histogrTitle),'_AverageCumulExcess']) = AverageCumulExcess;           
                currentTable.([char(currDataName),'_',char(histogrTitle),'_cumulExcess_hiCI']) = cumulExcess_hiCI;     
                currentTable.([char(currDataName),'_',char(histogrTitle),'_cumulExcess_loCI']) = cumulExcess_loCI;     
                currentTable.([char(currDataName),'_',char(histogrTitle),'_meanZoneWidth']) = meanZoneWidth;              
                currentTable.([char(currDataName),'_',char(histogrTitle),'_meanExcessDirection']) = meanExcessDirection;
 
                clear nelements centers cumulExcess zoneWidth ExcessDirection exist_Lumen exist_Adjacent 
                
            catch               
                currentTable.([char(currDataName),'_',char(histogrTitle),'_AverageCumulExcess']) = 0;           
                currentTable.([char(currDataName),'_',char(histogrTitle),'_cumulExcess_hiCI']) = 0;     
                currentTable.([char(currDataName),'_',char(histogrTitle),'_cumulExcess_loCI']) = 0;     
                currentTable.([char(currDataName),'_',char(histogrTitle),'_meanZoneWidth']) = -1;              
                currentTable.([char(currDataName),'_',char(histogrTitle),'_meanExcessDirection']) = 0;
                
                warning(['not found or error while processing: ',cnst.cachepath,currImage.tiffname,...
                '_distances-histogram-Data',char(currDataName),...
                'class_',strrep(char(histogrTitle),'_',' '),'.mat']);
            
            end
        end
    end
    
    
    %----------
    
    
    % save to master table
    if (first == 1)
     masterTable = struct2table(currentTable);
     first = 0;
    else 
     masterTable = [masterTable; struct2table(currentTable)];
    end
    

    else
       warning('omitted ROI because name is not == FullTumorValidated'); 
    end
    
     
end
end

writetable(masterTable,'QuantificationHotspotTable_FIRST_COHORT_2016-06-27.csv');
writetable(masterTable,'QuantificationHotspotTable_FIRST_COHORT_2016-06-27.xlsx');

