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

% ROTATE ROI
clear all, close all, format compact, clc

% -- add toolboxes to current path
addpath([pwd,'/subroutines'],'-end'); % my own subroutines
addpath([pwd,'/poly_stuff'],'-end');  % code from http://www.mathworks.com/matlabcentral/fileexchange/10391-fast-points-in-polygon-test
addpath([pwd,'/Yen'],'-end');  % code from http://www.mathworks.com/matlabcentral/fileexchange/10391-fast-points-in-polygon-test
addpath([pwd,'/registration'],'-end'); % registration code by CAW/JNK

cnst = getConstants(); % define constants

% additional constants. To do: move this block to getConstants()
cnst.Ki67folder = 'RegistrInput/Ki67/';
cnst.CD03folder = 'RegistrInput/CD03/';
cnst.cachefolder = '../cache_multivariate/';
cnst.thumbfolder = 'thumbs_multivariate/';
cnst.thumbViewFactor = 0.1;
cnst.regOutputFolder = 'RegistrOutput/';
% ------------------------

for iterateID = {'Smp006', 'Smp007'}
    

disp(['STARTING ',char(iterateID)]);
load([cnst.regOutputFolder,char(iterateID),'.mat']); % load dataset

% assume that a tumor region exists
TumorPoly_orig = IMGdataset.CD34.ROIpolygons{IMGdataset.CD34.TumorROInum};

% try: if a lumen region exists, extract the corresponding coordinates
try
    LumenPoly_orig = IMGdataset.CD34.ROIpolygons{IMGdataset.CD34.LumenROInum};
catch
    LumenPoly_orig = zeros(2,2);
end

% try if a adjacent region exists, extract the corresponding coordinates
try
    AdjacentPoly_orig = IMGdataset.CD34.ROIpolygons{IMGdataset.CD34.AdjacentROInum};
catch
    AdjacentPoly_orig = zeros(2,2);
end

%% PROCESS Ki67
% -----
% rotate ROI to fit Ki67 image / data
TumorPoly_Ki67 = TumorPoly_orig * IMGdataset.Ki67.ROIrot; % rotate
TumorPoly_Ki67(:,1) = TumorPoly_Ki67(:,1) + IMGdataset.Ki67.ROIxMove / cnst.thumbViewFactor; % move x
TumorPoly_Ki67(:,2) = TumorPoly_Ki67(:,2) + IMGdataset.Ki67.ROIyMove / cnst.thumbViewFactor; % move y

% ----
% rotate Ki67 brown coordinates to fit original ROI
IMGdataset.Ki67.brownCoordsRotated = IMGdataset.Ki67.brownCoords * IMGdataset.Ki67.ROIrotINV; % rotate all Ki67 pts back
IMGdataset.Ki67.brownCoordsRotated(:,1) = IMGdataset.Ki67.brownCoordsRotated(:,1) ...
    - IMGdataset.Ki67.ROIxMove / cnst.thumbViewFactor; % move x
IMGdataset.Ki67.brownCoordsRotated(:,2) = IMGdataset.Ki67.brownCoordsRotated(:,2) ...
    - IMGdataset.Ki67.ROIyMove / cnst.thumbViewFactor; % move y

% rotate Ki67 blue coordinates to fit original ROI
IMGdataset.Ki67.blueCoordsRotated = IMGdataset.Ki67.blueCoords * IMGdataset.Ki67.ROIrotINV; % rotate all Ki67 pts back
IMGdataset.Ki67.blueCoordsRotated(:,1) = IMGdataset.Ki67.blueCoordsRotated(:,1) ...
    - IMGdataset.Ki67.ROIxMove / cnst.thumbViewFactor; % move x
IMGdataset.Ki67.blueCoordsRotated(:,2) = IMGdataset.Ki67.blueCoordsRotated(:,2) ...
    - IMGdataset.Ki67.ROIyMove / cnst.thumbViewFactor; % move y

% -----
% extract Ki67 in-poly coords and perform statistics
[outStruct.Ki67.X,outStruct.Ki67.Y,outStruct.Ki67.probabilityMap,...
    outStruct.Ki67.densityOBS,outStruct.Ki67.sigLevelBonferroni] = ...
    hotspotProbability(IMGdataset.Ki67.brownCoordsRotated,...
    TumorPoly_orig,cnst);

% analyze distances for Ki67
quantify_CD03_Ki67_DistanceHistogram( ...
    IMGdataset.Ki67.brownCoordsRotated, ... % n-by-2-array of object coordinates
    [IMGdataset.Ki67.brownCoordsRotated; IMGdataset.Ki67.blueCoordsRotated], ... % n-by-2-array of all object coordinates
    imfinfo([cnst.img_folder,IMGdataset.CD34.fullfile]), ... % containing height and width of the original image
    IMGdataset.CD34.fullfile, ... % original image filename
    'Ki67_DATA', ... % type of coordinates used: e.g. Ki67
    IMGdataset.CD34.ROIpolygons, ... % contains ROI coordinates
    IMGdataset.CD34.ROInames, ... % contains ROI names
    cnst ) % contains constants

%% PROCESS CD03
% -----
% rotate ROI to fit CD3 image / data
TumorPoly_CD03 = TumorPoly_orig;
TumorPoly_CD03 = TumorPoly_CD03 * IMGdataset.CD03.ROIrot; % rotate
TumorPoly_CD03(:,1) = TumorPoly_CD03(:,1) + IMGdataset.CD03.ROIxMove / cnst.thumbViewFactor; % move x
TumorPoly_CD03(:,2) = TumorPoly_CD03(:,2) + IMGdataset.CD03.ROIyMove / cnst.thumbViewFactor; % move y

% ----
% rotate CD03 brown coordinates to fit original ROI
IMGdataset.CD03.brownCoordsRotated = IMGdataset.CD03.brownCoords * IMGdataset.CD03.ROIrotINV; % rotate all CD03 pts back
IMGdataset.CD03.brownCoordsRotated(:,1) = IMGdataset.CD03.brownCoordsRotated(:,1) ...
    - IMGdataset.CD03.ROIxMove / cnst.thumbViewFactor; % move x
IMGdataset.CD03.brownCoordsRotated(:,2) = IMGdataset.CD03.brownCoordsRotated(:,2) ...
    - IMGdataset.CD03.ROIyMove / cnst.thumbViewFactor; % move y

% rotate CD03 blue coordinates to fit original ROI
IMGdataset.CD03.blueCoordsRotated = IMGdataset.CD03.blueCoords * IMGdataset.CD03.ROIrotINV; % rotate all Ki67 pts back
IMGdataset.CD03.blueCoordsRotated(:,1) = IMGdataset.CD03.blueCoordsRotated(:,1) ...
    - IMGdataset.CD03.ROIxMove / cnst.thumbViewFactor; % move x
IMGdataset.CD03.blueCoordsRotated(:,2) = IMGdataset.CD03.blueCoordsRotated(:,2) ...
    - IMGdataset.CD03.ROIyMove / cnst.thumbViewFactor; % move y

% -----
% extract CD3 in-poly coords and perform statistics
[outStruct.CD03.X,outStruct.CD03.Y,outStruct.CD03.probabilityMap,...
    outStruct.CD03.densityOBS,outStruct.CD03.sigLevelBonferroni] = ...
    hotspotProbability(IMGdataset.CD03.brownCoordsRotated,...
    TumorPoly_orig,cnst);

% analyze distances for CD3
quantify_CD03_Ki67_DistanceHistogram( ...
    IMGdataset.CD03.brownCoordsRotated, ... % n-by-2-array of object coordinates
    [IMGdataset.CD03.brownCoordsRotated; IMGdataset.CD03.blueCoordsRotated], ... % n-by-2-array of all object coordinates
    imfinfo([cnst.img_folder,IMGdataset.CD34.fullfile]), ... % containing height and width of the original image
    IMGdataset.CD34.fullfile, ... % original image filename
    'CD03_DATA', ... % type of coordinates used: e.g. Ki67
    IMGdataset.CD34.ROIpolygons, ... % contains ROI coordinates
    IMGdataset.CD34.ROInames, ... % contains ROI names
    cnst ) % contains constants

%% PLOT RESULTS
% --------------------------------
% --------------------------------
% --- SHOW DENSITY AND HOTSPOT PROBABILITY OF POINTS IN POLYGONAL ROI


% --- PLOT CD03 DENSITY 
figure()
%imshow((CD03thumb(:,:,1:3))); 
hold on
contourf(((outStruct.CD03.X))*cnst.thumbViewFactor,...
    ((outStruct.CD03.Y))*cnst.thumbViewFactor,...
    outStruct.CD03.densityOBS,...
    cnst.contourlines,'LineStyle','none');
plot(TumorPoly_orig(:,1)*cnst.thumbViewFactor,TumorPoly_orig(:,2)*cnst.thumbViewFactor,'w-','LineWidth',2);
patch(LumenPoly_orig(:,1)*cnst.thumbViewFactor,LumenPoly_orig(:,2)*cnst.thumbViewFactor,...
    [166,189,219]/255,...
    'FaceAlpha',1,'EdgeColor','white','LineStyle','none');
patch(AdjacentPoly_orig(:,1)*cnst.thumbViewFactor,AdjacentPoly_orig(:,2)*cnst.thumbViewFactor,...
    [28,144,153]/255,...
    'FaceAlpha',1,'EdgeColor','white','LineStyle','none');
title('CD03 density','FontSize',8);
hold off
axis equal tight off;
colormap hot, colorbar
% decorations
suptitle(IMGdataset.ID);
set(gcf,'Color','w');
drawnow
print(gcf,'-dpng','-r600',[cnst.cachepath,IMGdataset.ID,...
    '_MULTI_01_CD03_DENSITY.png']); 


% --- PLOT CD03 HOTSPOTS
figure()
hold on
contourf(((outStruct.CD03.X))*cnst.thumbViewFactor,...
    ((outStruct.CD03.Y))*cnst.thumbViewFactor,...
    outStruct.CD03.probabilityMap,...
    cnst.contourlines,'LineStyle','none');
plot(TumorPoly_orig(:,1)*cnst.thumbViewFactor,TumorPoly_orig(:,2)*cnst.thumbViewFactor,'w-','LineWidth',2);
patch(LumenPoly_orig(:,1)*cnst.thumbViewFactor,LumenPoly_orig(:,2)*cnst.thumbViewFactor,...
    [166,189,219]/255,...
    'FaceAlpha',1,'EdgeColor','white','LineStyle','none');
patch(AdjacentPoly_orig(:,1)*cnst.thumbViewFactor,AdjacentPoly_orig(:,2)*cnst.thumbViewFactor,...
    [28,144,153]/255,...
    'FaceAlpha',1,'EdgeColor','white','LineStyle','none');
title('CD03 hotspot probability','FontSize',8);
hold off
axis equal tight off;
zValueSigLevel = norminv(1-outStruct.CD03.sigLevelBonferroni)
colormap(hot());
caxis([0 zValueSigLevel/0.67]); 
colorbar;
% decorations
suptitle(IMGdataset.ID);
set(gcf,'Color','w');
drawnow
print(gcf,'-dpng','-r600',[cnst.cachepath,IMGdataset.ID,...
    '_MULTI_02_CD03_HOTSPOT_Zsig',num2str(zValueSigLevel),'.png']); 
clear zValueSigLevel

% --- PLOT Ki67 DENSITY
figure()
%imshow((Ki67thumb(:,:,1:3))); 
hold on
contourf(((outStruct.Ki67.X))*cnst.thumbViewFactor,...
    ((outStruct.Ki67.Y))*cnst.thumbViewFactor,...
    outStruct.Ki67.densityOBS,...
    cnst.contourlines,'LineStyle','none');
plot(TumorPoly_orig(:,1)*cnst.thumbViewFactor,TumorPoly_orig(:,2)*cnst.thumbViewFactor,'w-','LineWidth',2);
patch(LumenPoly_orig(:,1)*cnst.thumbViewFactor,LumenPoly_orig(:,2)*cnst.thumbViewFactor,...
    [166,189,219]/255,...
    'FaceAlpha',1,'EdgeColor','white','LineStyle','none');
patch(AdjacentPoly_orig(:,1)*cnst.thumbViewFactor,AdjacentPoly_orig(:,2)*cnst.thumbViewFactor,...
    [28,144,153]/255,...
    'FaceAlpha',1,'EdgeColor','white','LineStyle','none');
title('Ki67 density','FontSize',8);
hold off
axis equal tight off;
colormap hot, colorbar
% decorations
suptitle(IMGdataset.ID);
set(gcf,'Color','w');
drawnow
print(gcf,'-dpng','-r600',[cnst.cachepath,IMGdataset.ID,...
    '_MULTI_03_Ki67_DENSITY.png']); 

% --- PLOT Ki67 HOTSPOTS
figure()
hold on
contourf(((outStruct.Ki67.X))*cnst.thumbViewFactor,...
    ((outStruct.Ki67.Y))*cnst.thumbViewFactor,...
    outStruct.Ki67.probabilityMap,...
    cnst.contourlines,'LineStyle','none');
plot(TumorPoly_orig(:,1)*cnst.thumbViewFactor,TumorPoly_orig(:,2)*cnst.thumbViewFactor,'w-','LineWidth',2);
patch(LumenPoly_orig(:,1)*cnst.thumbViewFactor,LumenPoly_orig(:,2)*cnst.thumbViewFactor,...
    [166,189,219]/255,...
    'FaceAlpha',1,'EdgeColor','white','LineStyle','none');
patch(AdjacentPoly_orig(:,1)*cnst.thumbViewFactor,AdjacentPoly_orig(:,2)*cnst.thumbViewFactor,...
    [28,144,153]/255,...
    'FaceAlpha',1,'EdgeColor','white','LineStyle','none');
title('Ki67 hotspot probability','FontSize',8);
hold off
axis equal tight off;
zValueSigLevel = norminv(1-outStruct.Ki67.sigLevelBonferroni)
colormap(hot());
caxis([0 zValueSigLevel/0.67]);
colorbar;
% decorations
suptitle(IMGdataset.ID);
set(gcf,'Color','w');
drawnow
print(gcf,'-dpng','-r600',[cnst.cachepath,IMGdataset.ID,...
    '_MULTI_04_Ki67_HOTSPOT_Zsig',num2str(zValueSigLevel),'.png']); 
clear zValueSigLevel

% --- PLOT CD34 DENSITY
figure()
%imshow((CD34thumb(:,:,1:3))); 
hold on;
contourf(((IMGdataset.CD34.X))*cnst.thumbViewFactor,...
    ((IMGdataset.CD34.Y))*cnst.thumbViewFactor,...
    IMGdataset.CD34.densityMap,...
    cnst.contourlines,'LineStyle','none');
title('CD34 density','FontSize',8);
hold on
plot(TumorPoly_orig(:,1)*cnst.thumbViewFactor,TumorPoly_orig(:,2)*cnst.thumbViewFactor,'w-','LineWidth',2);
patch(LumenPoly_orig(:,1)*cnst.thumbViewFactor,LumenPoly_orig(:,2)*cnst.thumbViewFactor,...
    [166,189,219]/255,...
    'FaceAlpha',1,'EdgeColor','white','LineStyle','none');
patch(AdjacentPoly_orig(:,1)*cnst.thumbViewFactor,AdjacentPoly_orig(:,2)*cnst.thumbViewFactor,...
    [28,144,153]/255,...
    'FaceAlpha',1,'EdgeColor','white','LineStyle','none');
hold off
axis equal tight off;
colormap hot, colorbar
% decorations
suptitle(IMGdataset.ID);
set(gcf,'Color','w');
drawnow
print(gcf,'-dpng','-r600',[cnst.cachepath,IMGdataset.ID,...
    '_MULTI_05_CD34_DENSITY.png']); 

% --- PLOT CD34 HOTSPOTS
figure()
hold on
contourf(((IMGdataset.CD34.X))*cnst.thumbViewFactor,...
    ((IMGdataset.CD34.Y))*cnst.thumbViewFactor,...
    IMGdataset.CD34.probabilityMap,...
    cnst.contourlines,'LineStyle','none');
plot(TumorPoly_orig(:,1)*cnst.thumbViewFactor,TumorPoly_orig(:,2)*cnst.thumbViewFactor,'w-','LineWidth',2);
patch(LumenPoly_orig(:,1)*cnst.thumbViewFactor,LumenPoly_orig(:,2)*cnst.thumbViewFactor,...
    [166,189,219]/255,...
    'FaceAlpha',1,'EdgeColor','white','LineStyle','none');
patch(AdjacentPoly_orig(:,1)*cnst.thumbViewFactor,AdjacentPoly_orig(:,2)*cnst.thumbViewFactor,...
    [28,144,153]/255,...
    'FaceAlpha',1,'EdgeColor','white','LineStyle','none');
title('CD34 hotspot probability','FontSize',8);
hold off
axis equal tight off;
zValueSigLevel = norminv(1-IMGdataset.CD34.sigLevelBonferroni)
colormap(hot());
caxis([0 zValueSigLevel/0.67]);
colorbar;
% decorations
suptitle(IMGdataset.ID);
set(gcf,'Color','w');
drawnow
print(gcf,'-dpng','-r600',[cnst.cachepath,IMGdataset.ID,...
    '_MULTI_06_CD34_HOTSPOT_Zsig',num2str(zValueSigLevel),'.png']); 
clear zValueSigLevel


end
