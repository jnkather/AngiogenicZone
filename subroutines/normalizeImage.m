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

function imageOut = normalizeImage(imageIn, opt)

    % convert image to double
	imageOut = double(imageIn);

    % process each channel
	for i=1:size(imageIn,3)

        % extract one channel at a time
		Ch = imageOut(:,:,i);

        % normalize output range to 0...1
		imageOut(:,:,i) = (imageOut(:,:,i)-min(Ch(:)))/(max(Ch(:)-min(Ch(:))));

        % optional: stretch histogram
        if strcmp(opt,'stretch')
        	imageOut(:,:,i) = imadjust(imageOut(:,:,i), ...
						stretchlim(imageOut(:,:,i)),[]);
						 % default option for stretchlim is [0.01,0.99]
        end
	end

end
