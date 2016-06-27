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

function [dataReg] = regThumbnails(picFix, picMov, VisOff, vargin)

%% precalculation etc.
if exist('VisOff') ==0
    VisOff = 'On';
end

picFix = picFix(:,:,1:3);
picMov = picMov (:,:,1:3);

switch VisOff % visualization part
    case 'Off'
        
    case 'On'
    close all
    figure(1),
    subplot(1,2,1), imshow(picFix)
    subplot(1,2,2), imshow(picMov)
end

%% make images to a mask
% two basic assumptions: i) a mask reduces the complexity to find features;
% ii) a mask is smaller an therefore reduce computaion efforts

[maskPicFix shrinkageFactor] = im2mask(picFix); 
[maskPicMov] = im2mask(picMov, shrinkageFactor); 

switch VisOff % visualization part
    case 'Off'
        
    case 'On'
    close all
    figure(1), imagesc(maskPicFix), colormap 'hot'
    figure(2), imagesc(maskPicMov), colormap 'hot'
end

%% application of feature based registration
% detection of the features
ptsOriginal  = detectSURFFeatures(maskPicFix);
ptsDistorted = detectSURFFeatures(maskPicMov);
[featuresOriginal,  validPtsOriginal]  = extractFeatures(maskPicFix,  ptsOriginal);
[featuresDistorted, validPtsDistorted] = extractFeatures(maskPicMov, ptsDistorted);
indexPairs = matchFeatures(featuresOriginal, featuresDistorted);
matchedOriginal  = validPtsOriginal(indexPairs(:,1));
matchedDistorted = validPtsDistorted(indexPairs(:,2));

switch VisOff % visualization part
    case 'Off'
        
    case 'On'
    figure;
    showMatchedFeatures(maskPicFix,maskPicMov,matchedOriginal,matchedDistorted);
        title('Putatively matched points (including outliers)');
end

%% application of feature based registration 
% generation of the tform
[tform, inlierDistorted, inlierOriginal] = estimateGeometricTransform(...
    matchedDistorted, matchedOriginal, 'similarity');

switch VisOff % visualization part
    case 'Off'
        
    case 'On'
    figure;
    showMatchedFeatures(maskPicFix,maskPicMov,inlierOriginal,inlierDistorted);
        title('Matching points (inliers only)');
            legend('ptsOriginal','ptsDistorted');
end

%% application of feature based registration / intensity based does not work
% adaption of the tform and application
tform.T(3,1)= tform.T(3,1).* (1/shrinkageFactor);
tform.T(3,2)= tform.T(3,2).* (1/shrinkageFactor);

outputView = imref2d(size(picFix(:,:,3)));
picReg= imwarp(picMov,tform,'OutputView',outputView);

%% prepartion of the results
dataReg.Thumb = picReg;
% x-translation
dataReg.xMove = tform.T(3,1);
% y-translation
dataReg.yMove = tform.T(3,2);
% z-rotation as angle
Tinv=tform.invert.T;
ss =Tinv(2,1);
sc = Tinv(1,1);
theta = atan2(ss, sc)*180 /pi;
dataReg.zRotation = theta;
dataReg.tform = tform;
dataReg.shrinkageFactor = shrinkageFactor;

%% visualization part

switch VisOff % visualization part
    case 'Off'
        
    case 'On'
    close all
    figure(1), imshowpair(rgb2gray(picFix),rgb2gray(picReg),'Scaling','joint'), ...
        title 'combination of both after registration'
    figure(2), imshowpair(rgb2gray(picFix),rgb2gray(picMov),'Scaling','joint'), ...
        title 'combination of both before registration'
    
end



end % function

