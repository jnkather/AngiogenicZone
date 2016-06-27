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
function [currentTable,h13] = getHotspotStatistics(vesselCoordsX,...
                                             vesselCoordsY,...
                                             RoiCoordsX,...
                                             RoiCoordsY,...
                                             DensityX,...
                                             DensityY,...
                                             DensityZ,...
                                             probabilityMap,...
                                             sigLevelBonferroni,...
                                             imgName,...
                                             RoiName,...
                                             ObjTable,cnst)

    % clean up density array and hotspot array
    cleanDensityZ = DensityZ(:);
    cleanDensityZ(cleanDensityZ==0) = [];
    
    cleanProbability = probabilityMap(:);
    cleanProbability(cleanProbability==0) = [];
    
    % return simple statistics
    currentTable.imgID = cellstr(imgName); % image name
    currentTable.RoiID = cellstr(RoiName); % ROI name
    currentTable.ovrllDensMean = mean(cleanDensityZ(:)); % density function global mean
    currentTable.ovrllDensSTD = std(cleanDensityZ(:));   % density function STD
    currentTable.ovrllHSProbMean = mean(cleanProbability(:)); % probability map global mean
    currentTable.ovrllHSProbSTD = std(cleanProbability(:));   % probability map STD 
    currentTable.zValueSigLevel = norminv(1-sigLevelBonferroni);

    % --- calculate advanced statistics
    % which Vessel Centroids are in the ROI (as defined by a polygon)? 
    vesselsIN = inpoly([vesselCoordsY(:),vesselCoordsX(:)],[RoiCoordsY(:),RoiCoordsX(:)]);
    X_vesselsIN = vesselCoordsX(vesselsIN);
    Y_vesselsIN = vesselCoordsY(vesselsIN);
    
    % save geometric data to table
    currentTable.meanVesselAreaROI = mean(ObjTable.Area(vesselsIN));
    currentTable.meanMajorAxisLengthROI = mean(ObjTable.MajorAxisLength(vesselsIN));
    currentTable.meanMinorAxisLengthROI = mean(ObjTable.MinorAxisLength(vesselsIN));
    currentTable.meanEccentricityROI = mean(ObjTable.Eccentricity(vesselsIN));
    currentTable.meanOrientationROI = mean(ObjTable.Orientation(vesselsIN));
    currentTable.meanEquivDiameterROI = mean(ObjTable.EquivDiameter(vesselsIN));
    currentTable.meanSolidityROI = mean(ObjTable.Solidity(vesselsIN));
    currentTable.meanPerimeterROI = mean(ObjTable.Perimeter(vesselsIN));
    
%     figure()
%     scatter(X_vesselsIN,Y_vesselsIN,12,'black','filled');
%     hold on
%     patch(RoiCoordsX,...
%           RoiCoordsY,'red',...
%           'FaceAlpha',0.1,'EdgeColor','black','LineWidth',2);    
      
    % which sampling points are in the ROI?
    gridIN = inpoly([DensityX(:),DensityY(:)],[RoiCoordsX,RoiCoordsY]);
    X_samplingPtsROI = DensityX(gridIN);
    Y_samplingPtsROI = DensityY(gridIN);
    Z_samplingPtsROI = probabilityMap(gridIN);
    
%     figure()
%     scatter(X_samplingPtsROI,Y_samplingPtsROI,1,'black','filled');
%     axis equal
    
    % which sampling points are part of significant hotspots?
    SIG = (probabilityMap>=(norminv(1-sigLevelBonferroni)));
    X_samplingPtsHS = DensityX(SIG); 
    Y_samplingPtsHS = DensityY(SIG); 
    Z_samplingPtsHS = probabilityMap(SIG);
    
%     h11=figure();
%     hold on
%     plot(DensityX(:),DensityY(:),'b.')
%     plot(X_samplingPtsROI,Y_samplingPtsROI,'c.');
%     plot(X_samplingPtsHS,Y_samplingPtsHS,'r.');
%     axis equal tight
%     %disp('Number of sampling points, number of points inside ROI');
%     title(['Hotspots and Object area for ',10,strrep(imgName,'_',' ')])
    
%     hold on
%     scatter(X_samplingPtsHS,Y_samplingPtsHS,3,'blue','filled');   
    
    % label hotspots in ascending order, from left to right, using
    % 4-connectivity
    BWImg = reshape(SIG,numel(unique(DensityX)),numel(unique(DensityY)));
    labelsHotspots = bwlabel(BWImg,4);
    SIG_labeled = reshape(labelsHotspots,[],1);
    HotspotCount = max(unique(SIG_labeled));
    
    %figure(); hold on
    for i=1:HotspotCount
        % trace boundaries
        [r,c] = ind2sub(size(labelsHotspots),find(labelsHotspots==i,1));
        sammlung{i}.contour = bwtraceboundary(BWImg,[r c],'W',8,Inf,'counterclockwise');
        %Plot the contour on the image.
        %plot(sammlung{i}.contour(:,2),sammlung{i}.contour(:,1),'g','LineWidth',5);
    end

    % create mask of boundary pixels
    HotspotBoundaries = BWImg & ~imerode( BWImg, [0 1 0; 1 1 1; 0 1 0] ); 

    % label boundaries
    LabeledBoundaries = HotspotBoundaries .* labelsHotspots;

%     figure ()
%     imshow(LabeledBoundaries,[]); 
%     colormap([0,0,0; lines()]);
%     title('Boundaries only');
    
    % remove background from SIG_labeled and labeled Boundaries
    LabeledBoundaries= LabeledBoundaries(:);
    LabeledBoundaries(SIG_labeled==0) = [];
    SIG_labeled(SIG_labeled==0) = [];
    
    % show the labeled hotspots and all the vessel centroids
    h13 = figure();
    hold on
     patch(RoiCoordsX,...
          RoiCoordsY,'red',...
          'FaceAlpha',0.1,'EdgeColor','black','LineWidth',2); 
    plot(X_vesselsIN,Y_vesselsIN,'.','Color',[.5 .5 .5]);
    scatter(X_samplingPtsHS,Y_samplingPtsHS,50,SIG_labeled,'filled')
    axis equal
    hold off
    title(['Labeled hotspots for ',10,strrep(imgName,'_',' '),...
        10,strrep(RoiName,'_',' ')])
    set(gcf,'Color','w');
    
    % plot the concave hull
    h14 = figure();
    axis equal
    hold on

    %scale coordinates
    scaleX = (max(DensityX(:))-min(DensityX(:)))/size(BWImg,1);
    scaleY = (max(DensityY(:))-min(DensityY(:)))/size(BWImg,2);
    XHalfstep = scaleX/2;
    YHalfstep = scaleY/2;
    numberVesselsInHotSpot = zeros(HotspotCount,1);

    for i=1:HotspotCount
        %create polygon

        x2 = min(DensityX(:))-XHalfstep+scaleX*sammlung{i}.contour(:,2);
        y2 = min(DensityY(:))-YHalfstep+scaleY*sammlung{i}.contour(:,1);

        VesselIsInHotSpot = inpolygon(X_vesselsIN,Y_vesselsIN,x2,y2);

        numberVesselsInHotSpot(i) = numel(X_vesselsIN(VesselIsInHotSpot));
     
        plot(X_samplingPtsHS(SIG_labeled==i),Y_samplingPtsHS(SIG_labeled==i),'c.','markersize',20);
        
        % save centers of hotspots
        HScenters(i,1) = mean(X_samplingPtsHS(SIG_labeled==i))
        HScenters(i,2) = mean(Y_samplingPtsHS(SIG_labeled==i))
 
        plot(X_vesselsIN(VesselIsInHotSpot),Y_vesselsIN(VesselIsInHotSpot),'r.','markersize',20);
        plot(x2,y2,'g-','LineWidth',5);
        
    end

    hold off
    
%     % plot hotspot centers -> to do: find 2nd order hotspots
%     figure()
%     plot(HScenters(:,1)*cnst.resizefactor,HScenters(:,2)*cnst.resizefactor,'r.')
%     hold on
%     plot(RoiCoordsX*cnst.resizefactor,RoiCoordsY*cnst.resizefactor,'b-');
%     axis equal tight
%     pause

    %--- calculate size etc.
    % Original scale = 0.495 microns per px

    % figure size in px
    maxX = max(DensityX(:))
    minX = min(DensityX(:))
    maxY = max(DensityY(:))
    minY = max(DensityY(:))
    rangeX = (max(DensityX(:))-min(DensityX(:)))
    rangeY = (max(DensityY(:))-min(DensityY(:)))
    rangeXsc = rangeX * 0.495 * 10^-3
    rangeYsc = rangeY * 0.495 * 10^-3
    
    sizeOfFigPx = (max(DensityX(:))-min(DensityX(:))) * (max(DensityY(:))-min(DensityY(:)));

    % convert fig size to square millimeters
    sizeOfFigMM2 = sizeOfFigPx * 0.495^2 * 10^-6;

    squareMMperHotspot=zeros(HotspotCount,1);
    for i=1:HotspotCount
        squareMMperHotspot(i) = (sum(SIG_labeled==i)/numel(DensityX(:)))*sizeOfFigMM2;
    end
    
    % write to summary struct
    currentTable.SamplingPtsROI = numel(gridIN(gridIN>0));
    currentTable.SamplingPtsAllHS2 = numel(SIG(SIG>0));
    currentTable.SamplingPtsTotal = numel(DensityX(:));
    currentTable.sizeOfFigPx = sizeOfFigPx;

    currentTable.sizeOfROI_MM2 = ( numel(gridIN(gridIN>0)) / numel(gridIN(:)) ) *sizeOfFigMM2;
    currentTable.sizeOfFig_MM2 = sizeOfFigMM2;
    currentTable.HotspotCount = HotspotCount;

    currentTable.mean_HS_Size = mean(squareMMperHotspot);
    currentTable.dev_HS_Size = std(squareMMperHotspot);

    vesselsPerMM = numberVesselsInHotSpot(:) ./ squareMMperHotspot(:);

    currentTable.overallMVDInROI = numel(vesselsIN(vesselsIN>0)) / currentTable.sizeOfROI_MM2;
    currentTable.meanMVDInHotspots = mean(vesselsPerMM);
    currentTable.stdMVDInHotspots  =std(vesselsPerMM);
    % currentTable.bandwidth  =bandwidth_CSRData(1);

    currentTable.HotspotsPerMM2 = HotspotCount / currentTable.sizeOfROI_MM2;
    
    % save advanced statistics to output struct
    currentTable.nHS = HotspotCount;
end

