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
function imageOut = SeparateStains(imageRGB, Matrix, varargin)

    % convert input image to double precision float
    imageRGB = double(imageRGB);

    % avoid log artifacts
    imageRGB = imageRGB + 2; % this is how scikit does it
    % imageRGB(imageRGB == 0) = 0.001; % this is how Fiji does it

    % perform color deconvolution: convert to OD then matrix multiplication
    imageOut = reshape(-log(imageRGB),[],3) * Matrix;
    imageOut = reshape(imageOut, size(imageRGB));     % reconstruct image

    % post-processing
    if nargin>2
        normalizeArgument = varargin{1};
    else
        normalizeArgument = 'stretch'; % default is 'stretch'
    end

    imageOut = normalizeImage(imageOut,normalizeArgument);    % normalize histogram
    imageOut = 1 - imageOut; % invert image
end
