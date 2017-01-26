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
function [xMean, loCI, hiCI, decis, pVal] = getMeanAndCI(dataContainer)

    % calculate confidence intervals
    xMean = mean(dataContainer);
    SEM = std(dataContainer)/sqrt(numel(dataContainer));    % standard error of the mean
    hiCI = xMean + 1.96 * SEM; % calculate upper 95% CI
    loCI = xMean - 1.96 * SEM; % calculate lower 95% CI
    
    % test the hypothesis that 0 comes from the observed distribution
    [decis, pVal] = ttest(dataContainer,0);
end

