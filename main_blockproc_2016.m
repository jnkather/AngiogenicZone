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
%% HEADER
% -- initialize program
clear all; close all; format compact; clc; 

% -- add toolboxes to current path
addpath([pwd,'/subroutines'],'-end'); % my own subroutines
addpath([pwd,'/poly_stuff'],'-end');  % code from http://www.mathworks.com/matlabcentral/fileexchange/10391-fast-points-in-polygon-test
addpath([pwd,'/Yen'],'-end');  % Yen thresholding

cnst = getConstants(); % define constants

% specify images (provide cell array)
imgArray = {'C2-Smp001-CD34.tif',...
            'C2-Smp002-CD34.tif',...
            'C2-Smp003-CD34.tif',...
            'C2-Smp004-CD34.tif',...
            'C2-Smp005-CD34.tif',...
            'C2-Smp006-CD34.tif',...
            'C2-Smp007-CD34.tif',...
            'C2-Smp008-CD34.tif',...
            'C2-Smp009-CD34.tif',...
            'C2-Smp010-CD34.tif',...
            'C2-Smp011-CD34.tif',...
            'C2-Smp012-CD34.tif',...
            'C2-Smp013-CD34.tif',...
            'C2-Smp014-CD34.tif',...
            'C2-Smp015-CD34.tif',...
            'C2-Smp016-CD34.tif',...
            'C2-Smp017-CD34.tif',...
            'C2-Smp018-CD34.tif',...
            'C2-Smp019-CD34.tif',...
            'C2-Smp020-CD34.tif',...
            'C2-Smp021-CD34.tif',...
            'C2-Smp022-CD34.tif',...
            'C2-Smp023-CD34.tif',...
            'C2-Smp024-CD34.tif',...
            'C2-Smp025-CD34.tif',...
            'C2-Smp026-CD34.tif',...
            'C2-Smp027-CD34.tif',...
            'C2-Smp028-CD34.tif',...
            'C2-Smp029-CD34.tif'};


% {'Smp001-CD34.tif',...
%             'Smp002-CD34.tif',...
%             'Smp003-CD34.tif',...
%             'Smp004-CD34.tif',...
%             'Smp005-CD34.tif',...
%             'Smp006-CD34.tif',...
%             'Smp007-CD34.tif',...
%             'Smp008-CD34.tif',...
%             'Smp009-CD34.tif',...
%             'Smp010-CD34.tif',...
%             'Smp011-CD34.tif',...
%             'Smp012-CD34.tif',...
%             'Smp013-CD34.tif',...
%             'Smp014-CD34.tif',...
%             'Smp015-CD34.tif',...
%             'Smp016-CD34.tif',...
%             'Smp017-CD34.tif',...
%             'Smp018-CD34.tif',...
%             'Smp019-CD34.tif',...
%             'Smp020-CD34.tif',...
%             'Smp021-CD34.tif',...
%             'Smp022-CD34.tif',...
%             'Smp023-CD34.tif',...
%             'Smp024-CD34.tif',...
%             'Smp025-CD34.tif',...
%             'Smp026-CD34.tif',...
%             'Smp027-CD34.tif',...
%             'Smp028-CD34.tif',...
%             'Smp029-CD34.tif',...
%             'Smp030-CD34.tif',...
%             'Smp031-CD34.tif',...
%             'Smp032-CD34.tif',...
%             'Smp033-CD34.tif',...
%             'Smp034-CD34.tif'};
        
        
%         fliplr({'Smp035-CD34.tif',...
%             'Smp036-CD34.tif',...
%             'Smp037-CD34.tif',...
%             'Smp038-CD34.tif',...
%             'Smp039-CD34.tif',...
%             'Smp040-CD34.tif',...
%             'Smp041-CD34.tif',...
%             'Smp043-CD34.tif',...
%             'Smp044-CD34.tif',...
%             'Smp045-CD34.tif',...
%             'Smp046-CD34.tif',...
%             'Smp047-CD34.tif',...
%             'Smp048-CD34.tif',...
%             'Smp049-CD34.tif',...
%             'Smp050-CD34.tif',...
%             'Smp051-CD34.tif',...
%             'Smp052-CD34.tif',...
%             'Smp053-CD34.tif',...
%             'Smp054-CD34.tif',...
%             'Smp056-CD34.tif',...
%             'Smp057-CD34.tif',...
%             'Smp058-CD34.tif',...
%             'Smp060-CD34.tif',...
%             'Smp061-CD34.tif',...
%             'Smp062-CD34.tif',...
%             'Smp063-CD34.tif',...
%             'Smp065-CD34.tif'});

        


% % caution: these images have to be 3-channel tiff images

% caution: MATLAB version:
% - whole slide analysis: use Matlab R2014a
% - distance maps: requires Matlab R2015a because of impixel and plotting

% --- Configuration block
% The analysis of whole slide images using MATLAB can take a while.
% Typically, the full procedure takes between 20 and 60 minutes per image.
% To make this analysis more practical, the different steps of the analysis
% can be performed independently of each other, provided that all raw data
% that is necessary for the current part is available. Here, logical
% variables can be used to define which block should be executed.
do.particleDetect = 0;  % analyze blood vessels in original TIFF image
do.requestROI = 0;      % ask the user to draw a ROI on each image
do.analysisInROI = 0;   % analyze the blood vessel distribution within each ROI
do.createThumb = 0;     % create a thumbnail image for each ROI
do.analyzeDistances = 1;  % analyze blood vessel distribution in the tumor

% -----------------------------------------------------------------------
if do.particleDetect
%% PART 1: Blood vessel detection
% iterate through all images for blood vessel particle detection
for imgNum = 1:numel(imgArray)
    
% specify image and read image metadata
currImage.tiffname = char(imgArray{imgNum});
currImage.metadata = imfinfo([cnst.path,currImage.tiffname]);  
disp(['vessel analysis | processing: ', currImage.tiffname]);   

try % try to load existing mask 
   load([cnst.cachepath,currImage.tiffname,'_vessel_measurements.mat']);
   warning(['Found measurements file.',...
       'Assuming this file contains a mask and image info.']);
catch % no mask found -> process image and create mask
    % perform block processing of the image to segment blood vessels
    fun = @(block_struct) rgb2vesselMask(block_struct.data,cnst);	
    currImage.Mask = blockproc([cnst.path,currImage.tiffname],...
                                       cnst.VesselBlckSize,fun,...
                                       'UseParallel',cnst.parallel,...
                                       'BorderSize',cnst.blockOverlap,...
                                       'PadPartialBlocks',true,...
                                       'PadMethod','symmetric');
                                   
    % cut padded areas. For multiframe Tiff files, Width and Heights are n-
    % dimensional arrays. In this case, select the first element. Observe 
    % that height and width are switched in this case.
    currImage.Mask = currImage.Mask(...
        1:(currImage.metadata(1).Height),...
        1:(currImage.metadata(1).Width));      
       
    % downscale mask
    currImage.MaskSmall = imresize(currImage.Mask, cnst.resizefactor);
    
    % perform morphological analysis of resized mask
    warning(['starting morphological analysis... this task is ',...
        'computationally intensive and may take a while']);
    tic, currImage.ObjTable = Mask2Table(currImage.MaskSmall,cnst); toc 
    
    % save measurement file
    save([cnst.cachepath,currImage.tiffname,'_vessel_measurements.mat'],...
                'currImage');
    imwrite(currImage.Mask,...
        [cnst.cachepath,currImage.tiffname,'_object_mask_control.jpeg'],...
        'Quality',60);
end
end % end iteration through all images
end
% -----------------------------------------------------------------------

% -----------------------------------------------------------------------
if do.requestROI
%% PART 2: Hotspot detection
clear currImage
% iterate through all images for ROI selection
for imgNum = 1:numel(imgArray)
close all; % close previous result windows
    
% specify image, read image metadata and load vessel measurements
currImage.tiffname = char(imgArray{imgNum});
currImage.metadata = imfinfo([cnst.path,currImage.tiffname]);  
load([cnst.cachepath,currImage.tiffname,'_vessel_measurements.mat']);
disp(['######',10,'hotspot analysis (if file exists, will overwrite) | processing: ', currImage.tiffname]);   

% get user input to define ROIs on current image
[currImage.polygonContainer, currImage.ROINameContainer] = ...
    getRoiSet(currImage, cnst);

if do.analysisInROI
% iterate through all ROIs for hotspot analysis
    for i=1:numel(currImage.ROINameContainer)
        if strcmp(char(currImage.ROINameContainer{i}),'FullTumorValidated')
            disp('Processing a FullTumorValidated ROI');
         
        % calculate hotspot probability
        [currMap.X,currMap.Y,currMap.probabilityMap,...
            currMap.densityOBS,currMap.sigLevelBonferroni] = ...
        hotspotProbability(currImage.ObjTable.Centroid/cnst.resizefactor,...
            currImage.polygonContainer{i},cnst);
        
        disp('finished calculating hotspot probability');
        
        if currMap.sigLevelBonferroni>0
            % show results and save results
            showHotspotResults(currImage,currMap, ...
                char(currImage.ROINameContainer{i}),...
                currImage.polygonContainer{i},cnst);
            disp('finished show hotspot results');
            % save hotspot measurement file for current ROI
            save([cnst.cachepath,currImage.tiffname,'_hotspot_measurements_',...
                char(currImage.ROINameContainer{i}),'.mat'], 'currMap');      
            disp('finished saving');
        else
            warning('did not save because sigLevel = 0');
        end
        else
            warning('did not save because ROI name ~= FullTumorValidated');
        end            
    end
end

end % end iteration through all images
end
% -----------------------------------------------------------------------

% -----------------------------------------------------------------------
if do.createThumb
%% PART 3: Save thumbnail
for imgNum = 1:numel(imgArray)
close all; % close previous result windows

% specify image, read image metadata and load vessel measurements
currImage.tiffname = char(imgArray{imgNum});
currImage.metadata = imfinfo([cnst.path,currImage.tiffname]);  
load([cnst.cachepath,currImage.tiffname,'_vessel_measurements.mat']);

disp(['##########',10,'processing: ',cnst.cachepath,...
    currImage.tiffname,'_vessel_measurements.mat']);
% find out if thumbnail exists
if(~numel(dir([cnst.cachepath,currImage.tiffname,'_Original.png'])))
    disp('file not found ... will create thumbnail');

% load ROIs
[currImage.polygonContainer, currImage.ROINameContainer] = ...
    getRoiSet(currImage, cnst);

figure(); 
warning('Loading image...');
tic
imshow([cnst.path,currImage.tiffname],...
    'Reduce',true,'InitialMagnification','fit');
toc
        
set(gcf,'Color','w');               % ...
hold on                             % plot all ROIs on top
for i=1:numel(currImage.polygonContainer)
    currPoly = currImage.polygonContainer{i};
    patch(currPoly(:,1),...
          currPoly(:,2),'red',...
          'FaceAlpha',0.1,'EdgeColor','black','LineWidth',2);          
end

% save the resulting figure as PNG
print(gcf,'-dpng','-r600',[cnst.cachepath,currImage.tiffname,...
    '_Original.png']);  
                              % finished this figure
else
    disp('thumbnail file found ... omitting');
end

end % end iteration through all images

end
% -----------------------------------------------------------------------

% -----------------------------------------------------------------------
if do.analyzeDistances
%% PART 5: Distribution of distances
clear currImage
% iterate through all images for distance analysis
for imgNum = 1:numel(imgArray)

% specify image, read image metadata and load vessel measurements
currImage.tiffname = char(imgArray{imgNum});
currImage.metadata = imfinfo([cnst.path,currImage.tiffname]);  
load([cnst.cachepath,currImage.tiffname,'_vessel_measurements.mat']);
disp(['######',10,'distance analysis | processing: ', currImage.tiffname]);   

% get user input to define ROIs on current image
[currImage.polygonContainer, currImage.ROINameContainer] = ...
    getRoiSet(currImage, cnst);

% --- start distance analysis
    for i=1:numel(currImage.ROINameContainer) % iterate through ROIs 
        % proceed only if the ROI is a 'FullTumorValidated' region
        if strcmp(char(currImage.ROINameContainer{i}),'FullTumorValidated')
        
        % prepare variables    
        disp('Processing a FullTumorValidated ROI');        
        Centroids = currImage.ObjTable.Centroid;
        currentROI =  currImage.polygonContainer{i}*cnst.resizefactor;
        
        % find in-polygon coordinates
        inP = inpolygon(Centroids(:,1),  Centroids(:,2), ...
                        currentROI(:,1), currentROI(:,2));             
        Centroids = Centroids(inP,:); % extract in-polygon coordinates
        
        % reset random number generator for reproducibility
        rng(cnst.randSeedDistances); 
        
        % create random point pattern
        tic
        disp('starting Monte Carlo experiment ...');
        for MCiterate = 1:cnst.MonteCarloIterations
            [CentroidsRANDX, CentroidsRANDY] = shuffleCoordsInPoly(Centroids(:,1), Centroids(:,2), ...
            currentROI(:,1),currentROI(:,2),  15 );
            Centroids_R{MCiterate} = [CentroidsRANDX, CentroidsRANDY]; % save the random points
        end
        toc 
        disp('finished Monte Carlo experiment ...');
        
        % preallocate
        currentAdjacentTissue = [];
        currentLumen = [];
            
        % iterate through ROINameContainer and find index of 
        % AdjacentTissueValidated and LumenValidated
        for k=1:numel(currImage.ROINameContainer)
            if strcmp(char(currImage.ROINameContainer{k}),'AdjacentTissueValidated')
                currentAdjacentTissue = currImage.polygonContainer{k}*cnst.resizefactor;
            elseif strcmp(char(currImage.ROINameContainer{k}),'LumenValidated')
                currentLumen = currImage.polygonContainer{k}*cnst.resizefactor;                
            end % end if 
        end % end for (look for AdjacentTissueValidated / LumenValidated)
        disp('identified ROIs');  

        %-- convert ROIs to binary masks
        % preallocate RAM for ROIs saved as binary masks
        Tissue_Masks = zeros(round(currImage.metadata(1).Height * cnst.resizefactor),...
                             round(currImage.metadata(1).Width * cnst.resizefactor),...
                             3);
                
        % process tumor region:
        % save the inverse of the tumor area as a binary mask
        Tissue_Masks(:,:,1) = ~poly2mask(currentROI(:,1),...
                currentROI(:,2),...
                round(currImage.metadata(1).Height * cnst.resizefactor),...
                round(currImage.metadata(1).Width * cnst.resizefactor));
        
        % process adjacent tissue region    
        try % try to save the adjacent tissue polygon as a binary mask
            Tissue_Masks(:,:,2) = poly2mask(currentAdjacentTissue(:,1),...
                        currentAdjacentTissue(:,2),...
                        round(currImage.metadata(1).Height * cnst.resizefactor),...
                        round(currImage.metadata(1).Width * cnst.resizefactor));
                    
             exist_Adjacent = 1;     
             disp('created binary mask for Adjacent');
        catch
            warning('NO ADJACENT TISSUE');
            exist_Adjacent = 0;
        end
        
        % process lumen region
        try % try to save the Lumen polygon as a binary mask
            Tissue_Masks(:,:,3) = poly2mask(currentLumen(:,1),...
                currentLumen(:,2),...
                round(currImage.metadata(1).Height * cnst.resizefactor),...
                round(currImage.metadata(1).Width * cnst.resizefactor));
            exist_Lumen = 1;
            disp('created binary mask for Lumen');
        catch
            warning('NO LUMEN');
            exist_Lumen = 0; 
        end
        
        % convert binary masks to distance matrix
        %% distance from class i: 
        % iterate through all classes and for each object
        % calculate the distance to each class
        distancesToClass = double(zeros(size(Centroids,1),3)); % preallocate
        distancesToClass_RAND = double(zeros(size(Centroids,1),cnst.MonteCarloIterations));  % preallocate
        for TumTisLum = 1:3 % iterate through FullTumor, AdjacentTissue, Lumen
            disp('starting iteration (TumTisLum), calculating distance transform');
            % create distance matrix for current class
            I_dist_class = bwdist(double(Tissue_Masks(:,:,TumTisLum)));
 
            % extract distances for observed data for current class
            temp_distances = impixel(I_dist_class,...
                        Centroids(:,1),Centroids(:,2));
            distancesToClass(:,TumTisLum) = double(temp_distances(:,1)) * cnst.distanceScalingMM;
            clear temp_distances
            
            % extract distances for expected (random) data for current class
            
            % preallocate
            distancesToClass_RAND = double(zeros(size(Centroids,1),...
                cnst.MonteCarloIterations));
            
            % iterate through all experiments
            disp('reading data of Monte Carlo experiments for current class...')
            tic
            for MCiterate = 1:cnst.MonteCarloIterations 
                % extract random coordinates of nth experiment
                currentRandCentroids = Centroids_R{MCiterate}; 
                temp_distances = impixel(I_dist_class,... % measure distances 
                        currentRandCentroids(:,1),currentRandCentroids(:,2)); 
                
                % distancesToClass_RAND is a MxN matrix:
                % M = number of points (given by size(Centroids,1))
                % N = number of experiments (given by cnst.MonteCarloIterations)
                distancesToClass_RAND(:,MCiterate) = double(temp_distances(:,1)) ...
                    * cnst.distanceScalingMM; % save distances of nth experiment
            end
            toc
            
            % save all N experiments to class container
            randDistancesContainer{TumTisLum} = distancesToClass_RAND;
            disp(['saved class results to container',10,'######']);
            
            clear temp_distances distancesToClass_RAND     
        end % end iteration through three classes of tissue types
        

        % display distance results
        for TumTisLum = 2:3 % iterate through class 2 (AdjTissue) and 
                               % class 3(Lumen) and plot distance to 
                               % AdjTissue and Lumen
                               
            % create the vector containing centers for the histogram                   
            xbinV = 0:cnst.distanceBinWidth:cnst.distanceMax;
            
            % create histogram data for observed data
            [nelements_OBS,centers] = hist(distancesToClass(:,TumTisLum),...
                            xbinV);
            
            % extract random distances for current class
            distancesToClass_RAND = randDistancesContainer{TumTisLum};
                                                
            % create histogram data for random data
            % iterate through all experiments
            for MCiterate = 1:cnst.MonteCarloIterations 
                [nelements_rand(:,MCiterate),centers] = ...
                    hist(distancesToClass_RAND(:,MCiterate), xbinV);
                % calculate excess vessels for current Monte Carlo experiment
                [ cumulExcessContainer(MCiterate), ...
                    zoneWidthContainer(MCiterate), ...
                    ExcessDirectionContainer(MCiterate), ~, ~ ] = ...
                    calculateExcessObjects( nelements_OBS, nelements_rand(:,MCiterate)', cnst);
            end

            % calculate global mean and confidence interval and other statistics    
            [AverageCumulExcess, cumulExcess_loCI, cumulExcess_hiCI, ttestDecis, ttestPVal] = ...
                getMeanAndCI(cumulExcessContainer);
            CumulExcessStd = std(cumulExcessContainer);
            meanZoneWidth = mean(zoneWidthContainer);
            meanExcessDirection = mean(ExcessDirectionContainer);
            
            clear distancesToClass_RAND % clean up
            
            % calculate mean of N Monte Carlo Experiments
            nelements_RAND_MEAN = mean(nelements_rand,2)';
            
            % show results of Monte Carlo distance distributions
            figure()
            subplot(1,2,1)
            plot(centers,nelements_rand,'k:')
            hold on; 
            plot(centers,nelements_RAND_MEAN,'b-','LineWidth',2)
            hold on;
            plot(centers,nelements_OBS,'r-','LineWidth',2);
            axTemp = axis();
            axTemp(1) = 0; axTemp(2) = cnst.distanceMax;
            axis(axTemp);
            xlabel('distance (mm)'); ylabel('aggregated blood vessel count');
                                   
            % calculate excess vessels for MEAN
           [ titleExcess, ~, ~, excess_1_MEAN, excess_2_MEAN ] = ...
                calculateExcessObjects( nelements_OBS, nelements_RAND_MEAN, cnst);
            
            % open up a subplot and show the object excess as a bar chart
            subplot(1,2,2)
            bar(centers,excess_1_MEAN,1); % plot the difference between OBS and RAND
            axis([0 cnst.distanceMax -300 300]);
            histogrTitle = char(cnst.ROIclassnames{TumTisLum});
            
            %legend({'excess vessels','baseline'});
            xlabel('distance (mm)'); ylabel('excess blood vessels');
              
            % plot decorations
            title([histogrTitle, 10, 'excess of the mean curve ', ...
                num2str(titleExcess),10,'mean of excess vectors ', num2str(AverageCumulExcess)]); % decorations
            % more plot decorations
            set(gcf,'Color','w');      
            suptitle([strrep(strrep(currImage.tiffname,'angiopredict','a-img'), '_', ' '),...
                10,'blood vessel distance from ...']);
            % save the resulting figure as PNG
            print(gcf,'-dpng','-r900',[cnst.cachepath,currImage.tiffname,...
                '_distances-normalized_HR_',histogrTitle,'.png']);
            
            % save results in cache
            save([cnst.cachepath,currImage.tiffname,...
            '_distances-histogram-class_',histogrTitle,'.mat'],...
            'nelements_OBS','nelements_RAND_MEAN','centers',...
            'meanZoneWidth','meanExcessDirection','AverageCumulExcess','CumulExcessStd',...
            'cumulExcess_loCI', 'cumulExcess_hiCI',...
            'ttestDecis', 'ttestPVal',...
            'exist_Lumen','exist_Adjacent'); % save in cache
        
            clear nelements_OBS nelements_RAND_MEAN centers meanZoneWidth AverageCumulExcess CumulExcessStd
            clear cumulExcess_loCI cumulExcess_hiCI ttestDecis ttestPVal meanExcessDirection titleExcess
        end
        
        disp('finished calculating distances');
        
        end % end processing FullTumorValidated       
        clear currentLumen currentAdjacentTissue currentROI 
    end % end iterating through ROInameContainer = end distance analysis
    clear currImage   
end % end iterating through images
end % end if do part 5
    
% -----------------------------------------------------------------------