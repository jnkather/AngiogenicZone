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
% getConstants: this function returns all constants for processing.
function cnst = getConstants()

    %% general settings
    cnst.img_folder = 'D:/Kather-LOCAL/My_clean_raw_WSI_data/2016-CRC_HLM_UMM_CD34/';  % all folder names end with /
    cnst.path = strcat(cnst.img_folder);
    cnst.cachepath = '../cache_multivariate/';
    cnst.resultsFolder = '../results_multivariate/'; 
    cnst.parallel = true;      % use parallel processing? default true
    cnst.resizefactor = 0.2;      % resize factor for saving masks
    cnst.DispOn = true;         % show output on screen? default true
    cnst.invert = true;         % invert image tile after processing? this may be neccessary
                                % depending on operating system and MATLAB
                                % version. default true
    cnst.saveResults = true;   % save results to files, default true
    cnst.noOverheat = false;     % pause after each segmentation run
    
    %% Blood vessel analysis
    cnst.VesselBlckSize = [5000 5000];  % tile size of original image
          % for TIFF images, tile size should be wisely chosen depending
          % on the TIFF tile size on the hard disk
    cnst.blockOverlap = [300 300];  % overlap of tiles in x & y dimension
    cnst.contourlines = 30;         % number of contour lines for plots
    cnst.sigLevel = 0.05; % level of significance before Bonferroni
    cnst.convType = 'Fiji';  % color deconv. type: 'Fiji' | 'optimal'
    cnst.seDilate = ones(9);    % smoothing kernel, default ones(9)
    cnst.seErode = ones(5);     % smoothing kernel, default ones(5)
    cnst.objThresh = 300;       % minimal object size in px, default 300
                                % object filtering occurs before resize
    cnst.RegionParams = {'Area','Centroid','EquivDiameter','Perimeter',...
         'Eccentricity','EulerNumber','Solidity','minorAxislength',...
         'majorAxislength','orientation'}; 
    % how many sampling points should the density function have? 
    cnst.DensityFunctionSampling = 1024; % has to be 2^N, default 1024
    cnst.KDE_kernel = [190,190];     % optimal bandwidth = 114/0.6 (ref.: our Oncotarget paper)
    cnst.saveMask = true;       % save mask as image
    cnst.normalizeArgu = 'stretch';
  
    %% analysis of distances
    cnst.ROIclassnames = {'NOT Tumor','Adjacent Tissue','Lumen'};
    cnst.distanceBinWidth = 0.2; % smoothing factor for distance histogram
    cnst.micronsPerPx = 0.5041; % microns per pixel in original image, First cohort: 0.495, Validation cohort: 0.5041
    cnst.distanceMax = 5; % max distance for histogram in mm, default 5
    cnst.distanceScalingMM = 1/cnst.resizefactor * cnst.micronsPerPx / 1000; % convert distance histogram x axis to mm
    cnst.randSeedDistances = 3; % for reproducibility of distance histograms, default 3
    cnst.MonteCarloIterations = 100; % Repeat the distance experiment N times, default 100
    
end