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
function [ X_OUT_CUMUL, Y_OUT_CUMUL ] = shuffleCoordsInPoly( X_IN, Y_IN, X_POLY, Y_POLY, n )
%shuffleCoordsInPoly shuffle 2D coordinates inside polygonal region
%   input: data set inside polygon, polygon coords
%   output: shuffled data set
%   n = number of iterations
%   this achieves csr within the ROI

    XMAX = max(X_POLY);
    YMAX = max(Y_POLY);
    XMIN = min(X_POLY);
    YMIN = min(Y_POLY);
    NUM = numel(X_IN);
    
    X_OUT_CUMUL = [1,1];
    Y_OUT_CUMUL = [1,1];
    
    for(i=1:n)
        X_OUT = (XMAX-XMIN) * rand(numel(X_IN),1) + XMIN;
        Y_OUT = (YMAX-YMIN) * rand(numel(Y_IN),1) + YMIN;

        IN = inpolygon(X_OUT,Y_OUT,X_POLY,Y_POLY);
        X_OUT = X_OUT(IN);
        Y_OUT = Y_OUT(IN);
        
        X_OUT_CUMUL = [X_OUT_CUMUL, X_OUT'];
        Y_OUT_CUMUL = [Y_OUT_CUMUL, Y_OUT'];
    end
    
    if(numel(X_OUT_CUMUL)<NUM)
        error('Shuffle failed');
    end    
    
    X_OUT_CUMUL = X_OUT_CUMUL(3:(NUM+2))';
    Y_OUT_CUMUL = Y_OUT_CUMUL(3:(NUM+2))';
    

end

