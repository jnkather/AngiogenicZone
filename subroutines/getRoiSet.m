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
% getRoiSet: this function tries to open an existing ROI Set for the current
% image. If a ROI Set is found, it is read and returned. If not, the user
% is prompted to define a ROI Set
function [polygonContainer, ROINameContainer] = getRoiSet(currImage,cnst)


    try     % try to load ROI set
        load([cnst.cachepath,currImage.tiffname,'_ROI_data.mat']);
        warning('A ROI set was found. Will proceed with this set and omit analysis.');
    
    catch   % if not found
        warning(['Did not find ',cnst.cachepath,currImage.tiffname,...
            '_ROI_data.mat','. Loading image...']);
        tic
        imshow([cnst.path,currImage.tiffname],...
            'Reduce',true,'InitialMagnification','fit');
        toc

        % draw a maximum of 10 polygons on the current image
        for(i=1:10)
            hCurrentPolygon = impoly();
            polygonContainer{i} = hCurrentPolygon.getPosition();
            wait(hCurrentPolygon)
            disp('Polygon was closed');

            % assign a name for the present ROI
            ROINameContainer{i} = inputdlg('Please assign a name to the selected ROI.');
                   
            % draw another polygon? if not, break loop
            hquest = questdlg('Do you want to draw another polygon?','User input required');
            if(strcmp(hquest,'No'))
               break
            elseif(strcmp(hquest,'Cancel'))
               break
            end
        end
        
        % save polygon container (coords and ROI names)
        save([cnst.cachepath,currImage.tiffname,'_ROI_data.mat'],...
            'polygonContainer','ROINameContainer');
    end

end