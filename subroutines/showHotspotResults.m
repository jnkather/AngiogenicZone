% created by JN Kather and CA Weis 2015-2016
% jakob.kather@nct-heidelberg.de
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
% showResults: this function creates plots and shows the results of the hotspot
% detection procedure
function showHotspotResults(currImage,currMap,currName,currentROI,cnst)

    % show objects as plot
    h1 = figure();
    scatter(currImage.ObjTable.Centroid(:,1),...
        currImage.ObjTable.Centroid(:,2),8,'filled');
    axis equal;
    title('Map of blood vessel centroids');
    set(gcf,'Color','w');

    % show density
    h2 = figure();
    contourf(currMap.X,currMap.Y,currMap.densityOBS,...
        cnst.contourlines,'LineStyle','none');
    hold on        
    % plot ROI on top
    patch(currentROI(:,1),currentROI(:,2),'red',...
           'FaceAlpha',0,'EdgeColor','white','LineWidth',2);        
    colormap(hot());
    colorbar;
    axis equal tight on;
    view(0,90);
    title(['blood vessel density map',10,'ROI ', ...
        currName]);
    set(gcf,'Color','w');
        
    % show hotspot probability map
    h3 = figure();
    % plot probability map as filled contour plot
    contourf(currMap.X,currMap.Y,currMap.probabilityMap,...
            cnst.contourlines,'LineStyle','none');
    hold on        
    % plot ROI on top
    patch(currentROI(:,1),currentROI(:,2),'red',...
           'FaceAlpha',0,'EdgeColor','white','LineWidth',2);         
    zValueSigLevel = norminv(1-currMap.sigLevelBonferroni)
    colormap(hot());
    caxis([0 zValueSigLevel/0.67]);
    colorbar;
    axis equal tight on;
    view(0,90);                    
    title(['angiogenic hotspot probability map',10,'ROI ', ...
        currName, ', level of significance at ',...
        num2str(zValueSigLevel), ' SD']);
    set(gcf,'Color','w');
    
   disp('starting to save...'); 
   if cnst.saveResults         % save the resulting figure as PNG
      print(h1,'-dpng','-r450',[cnst.resultsFolder,currImage.tiffname,...
           currName, '_01_SCATTER_', '.png']);  disp('done 1/3');
      print(h2,'-dpng','-r450',[cnst.resultsFolder,currImage.tiffname,...
           currName, '_02_DENSITY_', '.png']);  disp('done 2/3');
      print(h3,'-dpng','-r450',[cnst.resultsFolder,currImage.tiffname,...
           currName, '_02_HOTSPOTS_', '.png']); disp('done 3/3');
   end                         % finished this figure
        
end