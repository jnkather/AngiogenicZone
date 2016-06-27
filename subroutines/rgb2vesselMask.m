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
% processTile: this function is called for each tile of the image and
% performs all processing steps (thresholding, smoothing, filtering)
function imgOut = rgb2vesselMask(imgIn,cnst)
    
    tic
    disp('entering processTile');

    % if more than three channels, then remove additional channels
    if size(imgIn,3) > 3
        imgIn = imgIn(:,:,1:3);   
    end
    
    % get color deconvolution matrix for current tile
	[~,RGBtoHDAB] = getConversionMatrix(cnst.convType,imgIn);

    % separate stains = perform color deconvolution
    imageHDAB = SeparateStains(imgIn, RGBtoHDAB, cnst.normalizeArgu);
    
    % use difference of H and DAB channel for further processing
    imgDAB = imageHDAB(:,:,2);       
    imgH = imageHDAB(:,:,1);  
    
    % perform Yen segmentation of DAB
    numLevels = max(unique(imgDAB(:)));
    histogram = imhist(imgDAB,numLevels); 
    level = double(Yen(histogram))/255;
    imgMaskD = ~im2bw(imgDAB,level);

    % perform Yen segmentation of H
    numLevelsH = max(unique(imgH(:)));
    histogramH = imhist(imgH,numLevelsH); 
    levelH = double(Yen(histogramH))/255;
    imgMaskH = ~im2bw(imgH,levelH);
        
    % combine masks: remove blue objects
    imgMaskH = imerode(imdilate(imgMaskH,cnst.seDilate),cnst.seErode);
    imgMask = imgMaskD & ~imgMaskH;
  
    % dilate and erode  
    imgMask = imerode(imdilate(imgMask,cnst.seDilate),cnst.seErode);

    % remove particles smaller than objThresh
    imgOut = bwareaopen(imgMask, cnst.objThresh);

    % optional: pause to avoid core overheat
    elapsedTime = toc;
    if cnst.noOverheat, waitFor(elapsedTime,10); end
    
    disp('Completed blood vessel segmentation of current tile.');    
end