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
function [ cumulExcess, zoneWidth, ExcessDirection, excess_1, excess_2 ] = ...
    calculateExcessObjects( nelements_OBS, nelements_RAND, cnst)

    % prepare two types of object excess vectors            
    excess_1 = nelements_OBS(1:end) - nelements_RAND(1:end); % object excess starting at 1
    excess_2 = nelements_OBS(2:end) - nelements_RAND(2:end); % object excess starting at 2

     % try to extract the area and width of the first peak in the
            % histogram of the object excess
            try
                % disregard the first element
                excNEG = find(excess_2 < 0); % find the first negative 
                excPOS = find(excess_2 > 0); % find the first positive

                % calculate excess
                if excNEG(1) < excPOS(1)
                    % NEGATIVE EXCESS
                    cumulExcess = sum(excess_1(1:(excPOS(1))));
                    zoneWidth = (excPOS(1)) * cnst.distanceBinWidth;
                    ExcessDirection = -1;               
                elseif excNEG(1) > excPOS(1)
                    % POSITIVE EXCESS
                    cumulExcess = sum(excess_1(1:(excNEG(1))));
                    zoneWidth = (excNEG(1)) * cnst.distanceBinWidth;
                    ExcessDirection = 1;
                end
                
            catch
                cumulExcess = 0;
                zoneWidth = 0;
                ExcessDirection = 0;
                warning('histogram analysis failed ... ROI missing?');
            end

end

