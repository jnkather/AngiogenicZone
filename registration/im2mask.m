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
function [picOut, shrinkageFactor, vargout] = im2mask(picIn, shrinkageFactor, VisOn, vargin)


%% precalculations / perparations
if exist('VisOn') ==0
    VisOn = 'Off';
end

%%
picGray =256-rgb2gray(picIn);
picBw = im2bw(picGray, 0.10);

%% morphological operations
strelObject = strel('disk', 3); picBw = imdilate(picBw, strelObject);
strelObject = strel('disk', 5); picBw = imerode(picBw, strelObject);
strelObject = strel('disk', 10); picBw = imclose(picBw, strelObject);
strelObject = strel('disk', 15); picBw = imopen(picBw, strelObject);

%figure(1), imshow(picBw)
%figure(2), imshow(picMov)

%% prepare the output 

% loop to reduce the image size / computation time reduction
% this loop can calculate its own shrinkageFactor or use an input
[sizePicBwY sizePicBwX sizePicBwZ] = size(picBw);

if exist('shrinkageFactor') ==0 % loop to calculate
    
    if sizePicBwY > sizePicBwX 
        shrinkageFactor = 500 / sizePicBwY;
    else
        shrinkageFactor = 500 / sizePicBwX;
    end
    picBwSmall = imresize (picBw, shrinkageFactor);
    
else % loop to just use it
     picBwSmall = imresize (picBw, shrinkageFactor);
end
%% distance measurement applied
picDist = bwdist(~picBwSmall, 'euclidean');

switch VisOn
    
    case 'On'
        imagesc(picDist), colormap 'hot', colorbar
    case 'Off'
end

%% prepare results
picOut = double(picBwSmall);

end

