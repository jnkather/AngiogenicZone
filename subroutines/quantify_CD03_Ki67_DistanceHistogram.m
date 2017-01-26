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
function quantifyDistanceHistogram( ...
    Centroids_OBS, ... % n-by-2-array of object coordinates
    Centroids_CTRL, ... % n-by-2-array of control coordinates
    metadata, ... % containing height and width of the original image
    tiffname, ... % original image filename
    currDataName, ... % type of coordinates used: e.g. Ki67
    polygonContainer, ... % contains ROI coordinates
    ROINameContainer, ... % contains ROI names
    cnst ) % contains constants
    
% --- start distance analysis
    for i=1:numel(ROINameContainer) % iterate through ROIs 
        % proceed only if the ROI is a 'FullTumorValidated' region
        if strcmp(char(ROINameContainer{i}),'FullTumorValidated')
        % ------------------
        
        % prepare variables    
        disp('Processing a FullTumorValidated ROI');        
        currentROI =  polygonContainer{i}*cnst.resizefactor;
        Centroids_OBS = Centroids_OBS *cnst.resizefactor;
        Centroids_CTRL = Centroids_CTRL *cnst.resizefactor;
        % create the vector containing centers for the histogram                   
        xbinV = 0:cnst.distanceBinWidth:cnst.distanceMax;
        
        % find in-polygon coordinates
        inP_OBS = inpolygon(Centroids_OBS(:,1),  Centroids_OBS(:,2), ...
                        currentROI(:,1), currentROI(:,2));                       
        Centroids_OBS = Centroids_OBS(inP_OBS,:); % extract in-polygon coordinates
        
        inP_CTRL = inpolygon(Centroids_CTRL(:,1),  Centroids_CTRL(:,2), ...
                        currentROI(:,1), currentROI(:,2));   
        Centroids_CTRL = Centroids_CTRL(inP_CTRL,:);
           
        % preallocate
        currentAdjacentTissue = [];
        currentLumen = [];
        
        % iterate through ROINameContainer and find index of 
        % AdjacentTissueValidated and LumenValidated
        for k=1:numel(ROINameContainer)
            if strcmp(char(ROINameContainer{k}),'AdjacentTissueValidated')
                currentAdjacentTissue = polygonContainer{k}*cnst.resizefactor;
            elseif strcmp(char(ROINameContainer{k}),'LumenValidated')
                currentLumen = polygonContainer{k}*cnst.resizefactor;                
            end % end if 
        end % end for (look for AdjacentTissueValidated / LumenValidated)
        
        %-- convert ROIs to binary masks
        % preallocate RAM for ROIs saved as binary masks
        Tissue_Masks = zeros(round(metadata(1).Height * cnst.resizefactor),...
                             round(metadata(1).Width * cnst.resizefactor), 3);
                
        % process tumor region:
        % save the inverse of the tumor area as a binary mask
        Tissue_Masks(:,:,1) = ~poly2mask(currentROI(:,1),...
                currentROI(:,2),...
                round(metadata(1).Height * cnst.resizefactor),...
                round(metadata(1).Width * cnst.resizefactor));
        
        % process adjacent tissue region    
        try % try to save the adjacent tissue polygon as a binary mask
            Tissue_Masks(:,:,2) = poly2mask(currentAdjacentTissue(:,1),...
                        currentAdjacentTissue(:,2),...
                        round(metadata(1).Height * cnst.resizefactor),...
                        round(metadata(1).Width * cnst.resizefactor));    
             exist_Adjacent = 1;       
        catch
            warning('NO ADJACENT TISSUE');
            exist_Adjacent = 0;
        end
        
        % process lumen region
        try % try to save the Lumen polygon as a binary mask
            Tissue_Masks(:,:,3) = poly2mask(currentLumen(:,1),...
                        currentLumen(:,2),...
                        round(metadata(1).Height * cnst.resizefactor),...
                        round(metadata(1).Width * cnst.resizefactor));
                    exist_Lumen = 1;
        catch
            warning('NO LUMEN');
            exist_Lumen = 0; 
        end
        
        %% STATISTICS
        
        % perform the Monte Carlo experiment
        disp('starting Monte Carlo experiment');
        tic
        % reset random number generator for reproducibility
        rng(cnst.randSeedDistances);  
        % create N independent sets of random points
        % this is different than the random point generation for blood
        % vessels because here, the random points are constrained to the
        % locations of cells (blue or brown).        
        for MCiterate = 1:cnst.MonteCarloIterations 
            R_idx = randperm(size(Centroids_CTRL,1),size(Centroids_OBS,1));
            Centroids_R{MCiterate} = Centroids_CTRL(R_idx,:); % randomly draw k points
        end
        toc
        disp('created random coordinates');
               
        %% distance from class i: 
        % iterate through all classes and for each object
        % calculate the distance to each class
        distancesToClass = double(zeros(size(Centroids_OBS,1),3)); % preallocate
        distancesToClass_RAND = double(zeros(size(Centroids_OBS,1),cnst.MonteCarloIterations));  % preallocate
        for TumTisLum = 2:3 % iterate through FullTumor, AdjacentTissue, Lumen
            disp('starting iteration (TumTisLum), calculating distance transform');
            % create distance matrix for current class
            I_dist_class = bwdist(double(Tissue_Masks(:,:,TumTisLum)));
 
            % extract distances for observed data for current class
            temp_distances = impixel(I_dist_class,...
                        Centroids_OBS(:,1),Centroids_OBS(:,2));
            distancesToClass(:,TumTisLum) = double(temp_distances(:,1)) * cnst.distanceScalingMM;
            clear temp_distances
            
            % extract distances for expected (random) data for current class
            
            % preallocate
            distancesToClass_RAND = double(zeros(size(Centroids_OBS,1),...
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
            plot(centers,nelements_rand,'k-')
            hold on; 
            plot(centers,nelements_RAND_MEAN,'b-','LineWidth',2)
            hold on;
            plot(centers,nelements_OBS,'r-','LineWidth',2);
            axTemp = axis();
            axTemp(1) = 0; axTemp(2) = cnst.distanceMax;
            axis(axTemp);
            xlabel('distance (mm)'); ylabel('aggregated count');
                                   
            % calculate excess vessels for MEAN
           [ titleExcess, ~, ~, excess_1_MEAN, excess_2_MEAN ] = ...
                calculateExcessObjects( nelements_OBS, nelements_RAND_MEAN, cnst);
            
            % open up a subplot and show the object excess as a bar chart
            subplot(1,2,2)
            bar(centers,excess_1_MEAN,1); % plot the difference between OBS and RAND
            axis([0 cnst.distanceMax -5000 5000]);
            histogrTitle = char(cnst.ROIclassnames{TumTisLum});
            
            %legend({'excess vessels','baseline'});
            xlabel('distance (mm)'); ylabel(['excess ', currDataName]);
              
            % plot decorations
            title([histogrTitle, 10, 'excess of the mean curve ', ...
                num2str(titleExcess),10,'mean of excess vectors ', num2str(AverageCumulExcess)]); % decorations
            % more plot decorations
            set(gcf,'Color','w');      
            suptitle([strrep(strrep(tiffname,'angiopredict','a-img'), '_', ' '),...
                10,'object distance from ...']);
            % save the resulting figure as PNG
            print(gcf,'-dpng','-r900',[cnst.cachepath,tiffname,...
                currDataName,'-ROI_', histogrTitle,'_distances-normalized_HR_.png']);
            
            % save results in cache
            save([cnst.cachepath,tiffname,...
            '_distances-histogram-Data',currDataName,'class_',histogrTitle,'.mat'],...
            'AverageCumulExcess','cumulExcess_loCI','cumulExcess_hiCI','CumulExcessStd','meanZoneWidth','meanExcessDirection','ttestPVal',...
            'exist_Lumen','exist_Adjacent'); % save in cache
        
            clear nelements_OBS nelements_RAND_MEAN centers meanZoneWidth AverageCumulExcess CumulExcessStd
            clear cumulExcess_loCI cumulExcess_hiCI ttestDecis ttestPVal meanExcessDirection titleExcess
        end
        
        disp('finished calculating distances');        
        
        
        % ------------------
        end % end processing FullTumorValidated       
        clear currentLumen currentAdjacentTissue currentROI 
    end % end iterating through ROInameContainer = end distance analysis
end

