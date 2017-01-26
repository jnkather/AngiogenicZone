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
% hotspotProbability: this function uses the dataset of all blood vessel
% centroids as input and returns a hotspot probability map according to the
% algorithm from Kather et al., Oncotarget 2015, PMID 26061817

function [X,Y,probabilityMap,densityOBS,sigLevelBonferroni] = ...
    hotspotProbability(Centroids,currentROI,cnst)
    
    % extract in-polygon coordinates
    inP = inpolygon(Centroids(:,1),  Centroids(:,2), ...
                    currentROI(:,1), currentROI(:,2));
    
    if sum(inP) > 0
    
    disp([num2str(sum(inP)), ' sampling points out of a total of ', num2str(numel(inP)),...
        ' sampling points ',...
        ' are inside the current ROI']);
    
    % remove all centroids outside ROI
    Centroids(~inP,:) = [];
     
    % read min and max value
    MIN_XY = [min(currentROI(:,1)),min(currentROI(:,2))];
    MAX_XY = [max(currentROI(:,1)),max(currentROI(:,2))];
    BandwidthScaling = MAX_XY-MIN_XY;
    
    t_override = ([mean(cnst.KDE_kernel),mean(cnst.KDE_kernel)]./(BandwidthScaling)).^2;
    
    %% DENSITY OF OBSERVERD PATTERN VS. CSR
    % KDE of original pattern
    [bandwidthOBS,densityOBS,X,Y]=kde2d_JNK(Centroids,t_override,...
        cnst.DensityFunctionSampling,MIN_XY,MAX_XY);
    
    % create random pattern
     [Centroids_CSR_X,Centroids_CSR_Y] = shuffleCoordsInPoly( Centroids(:,1), Centroids(:,2), ...
        currentROI(:,1),currentROI(:,2),  15 );
    
    Centroids_CSR = [Centroids_CSR_X,Centroids_CSR_Y];

    % KDE of random pattern
    [bandwidthCSR,densityCSR,X2,Y2]=kde2d_JNK(Centroids_CSR,t_override,...
        cnst.DensityFunctionSampling,MIN_XY,MAX_XY);

    % only consider density function within the ROI
    % extract KDE grid points in ROI
    gridIN = inpolygon(X(:),Y(:),currentROI(:,1),currentROI(:,2));
    
%     % CALCULATE LEVEL OF SIGNIFICANCE
    sigLevelBonferroni = cnst.sigLevel / sum(gridIN);
    
    disp(['calculating hotspot probability... density function has ',...
        num2str(numel(densityOBS)), ' elements']);
    
    % normalize and return density
    probabilityMap = (densityOBS - mean(densityCSR(gridIN))) / ...
                            std(densityCSR(gridIN));
    probabilityMap(probabilityMap<0)=0;      
    
    else
        warning('There are no objects within the current ROI... omitting');
        X = 0;
        Y = 0;
        probabilityMap = 0;
        densityOBS = 0; 
        sigLevelBonferroni = 0;
    end

end

