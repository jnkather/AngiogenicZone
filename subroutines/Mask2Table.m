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
% Mask2Table: this function uses a full image object mask
% as input and returns the region properties of all objects in this 
% mask as a table
function ObjTable = Mask2Table(MaskIn,cnst)
    % extract object properties, return summary table
    Cstruct = bwconncomp(MaskIn,8); 
    disp('completed bwconncomp... startiong regionprops');
    ObjectInfo = regionprops(Cstruct, cnst.RegionParams);      
    ObjTable = struct2table(ObjectInfo);
end