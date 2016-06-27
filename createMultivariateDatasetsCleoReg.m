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
function createMultivariateDatasetsCleoReg()

    addpath([pwd,'/registration'],'-end'); % registration code by CAW/JNK

    % define constants
    cnst.Ki67folder = 'RegistrInput/Ki67/';
    cnst.CD03folder = 'RegistrInput/CD03/';
    cnst.cachefolder = '../cache_multivariate/';
    cnst.thumbfolder = 'thumbs_multivariate/';
    cnst.thumbViewFactor = 0.1;
    cnst.regOutputFolder = 'RegistrOutput/';
    
    % create dataset files - call the function for each image
    
    createDatasetFile('Smp026', ... % image ID
                      'Smp026-CD34.tif', ... % CD34 file
                      'Smp026-CD34.tif-THUMB10.tif', ... % CD34 thumbnail file
                      685, 763, cnst); % OMERO IDs for Ki67 and CD03 files

end


% --------------------------------------------------
% --------------------------------------------------
% --------------------------------------------------

%% master function
function createDatasetFile(datasetID_STR, CD34_tiffname_STR, CD34_thumbname_STR, ...
                            OMERO_ID_Ki67_INT, OMERO_ID_CD03_INT, cnst)
                        
%     try                    
    % --- start new dataset
    IMGdataset.ID = datasetID_STR;
    
    % pointers to full image files
    IMGdataset.CD34.fullfile = CD34_tiffname_STR;
    IMGdataset.Ki67.fullfile = '';
    IMGdataset.CD03.fullfile = '';
    
    % pointers to thumbnail image files
    IMGdataset.CD34.thumbfile = [cnst.thumbfolder,...
        CD34_thumbname_STR];
   
    % pointers to measurement files
    IMGdataset.CD03.measurementFile = [cnst.CD03folder,'Slide',num2str(OMERO_ID_CD03_INT),'.mat'];
    IMGdataset.Ki67.measurementFile = [cnst.Ki67folder,'Slide',num2str(OMERO_ID_Ki67_INT),'.mat'];
        
    % load CD34, CD03 and Ki67 and save result
    IMGdataset = loadSavedData(IMGdataset, cnst);
    
    % align thumbnails and save rotation and translation results
    IMGdataset = CalcROIOffset(IMGdataset,cnst);
    
    % write to file
    saveDataset(IMGdataset,cnst); 
    % end dataset
%     catch
%         warning(['ERROR IN ',datasetID_STR]);
%     end
end

%% LOAD FUNCTION ----

function IMGdataset = loadSavedData(IMGdataset,cnst)
    
    % create predetermined file names
    IMGdataset.CD03.thumbfile = [cnst.thumbfolder,IMGdataset.ID,'-CD03-THUMB10.tif'];
    IMGdataset.Ki67.thumbfile = [cnst.thumbfolder,IMGdataset.ID,'-Ki67-THUMB10.tif'];
    
    IMGdataset.CD34.measurementFileVES = [cnst.cachefolder, ... 
        IMGdataset.CD34.fullfile,'_vessel_measurements.mat'];
    IMGdataset.CD34.measurementFileHS = [cnst.cachefolder, ... 
        IMGdataset.CD34.fullfile,'_hotspot_measurements_FullTumorValidated.mat'];
    IMGdataset.CD34.measurementFileROI =[cnst.cachefolder, ... 
        IMGdataset.CD34.fullfile,'_ROI_data.mat'];
    
    % --------------------------------------------------------
    % --- load CD34 hotspot data
    load(IMGdataset.CD34.measurementFileHS)
    
    % save map to IMGdataset structure
    IMGdataset.CD34.X = currMap.X;
    IMGdataset.CD34.Y = currMap.Y;
    IMGdataset.CD34.probabilityMap = currMap.probabilityMap;
    IMGdataset.CD34.densityMap = currMap.densityOBS;
    IMGdataset.CD34.sigLevelBonferroni = currMap.sigLevelBonferroni;

    % --- load CD34 object data
    load(IMGdataset.CD34.measurementFileVES)
    
    % save blood vessel coordinates to IMGdataset structure   
    IMGdataset.CD34.VesselCoords = currImage.ObjTable.Centroid;
    
    % --------------------------------------------------------
    % --- load CD34 ROI data
    load(IMGdataset.CD34.measurementFileROI)
    IMGdataset.CD34.ROIpolygons = polygonContainer;
    IMGdataset.CD34.ROInames = ROINameContainer;
    
    % find FullTumorValidated
    for ROIi = 1:numel(IMGdataset.CD34.ROInames)
        if strcmp(IMGdataset.CD34.ROInames{ROIi},'FullTumorValidated')
            IMGdataset.CD34.TumorROInum = ROIi;
        end
    end 
    
    % find LumenValidated
    for ROIi = 1:numel(IMGdataset.CD34.ROInames)
        if strcmp(IMGdataset.CD34.ROInames{ROIi},'LumenValidated')
            IMGdataset.CD34.LumenROInum = ROIi;
        end
    end 
 
    % find AdjacentTissueValidated
    for ROIi = 1:numel(IMGdataset.CD34.ROInames)
        if strcmp(IMGdataset.CD34.ROInames{ROIi},'AdjacentTissueValidated')
            IMGdataset.CD34.AdjacentROInum = ROIi;
        end
    end 
    
    % stop only if no Tumor ROI found
    if ~IMGdataset.CD34.TumorROInum
        error('no FullTumorValidated found');
    end
    
    % --------------------------------------------------------
    % --- load CD3 data
    load(IMGdataset.CD03.measurementFile)
    IMGdataset.CD03.blueCoords = DataOut.blue;
    IMGdataset.CD03.brownCoords = DataOut.brown;

    % --------------------------------------------------------
    % --- load Ki67 data
    load(IMGdataset.Ki67.measurementFile)
    IMGdataset.Ki67.blueCoords = DataOut.blue;
    IMGdataset.Ki67.brownCoords = DataOut.brown;
    
end

%% SAVE FUNCTION ----

function saveDataset(IMGdataset,cnst)
    save([cnst.regOutputFolder,IMGdataset.ID,'.mat'],'IMGdataset');
    disp(['saved ',IMGdataset.ID]);
end


%% get thumbnails and rotate them

function IMGdataset = CalcROIOffset(IMGdataset,cnst)

% load CD34 as fix image
picFix = imread(IMGdataset.CD34.thumbfile);

% register CD34 and Ki67
picMov = imread(IMGdataset.Ki67.thumbfile);
[DataReg] = regThumbnails(picFix, picMov, 'Off');

IMGdataset.Ki67.xMove = DataReg.xMove;
    IMGdataset.Ki67.ROIxMove = - IMGdataset.Ki67.xMove;
IMGdataset.Ki67.yMove = DataReg.yMove;
    IMGdataset.Ki67.ROIyMove = - IMGdataset.Ki67.yMove;
IMGdataset.Ki67.zRotation = DataReg.zRotation;

IMGdataset.Ki67.thumfileReg = [cnst.regOutputFolder, IMGdataset.ID,'-Ki67-ThumbReg.tif' ];
imwrite(DataReg.Thumb, IMGdataset.Ki67.thumfileReg);

IMGdataset.Ki67.ROIrot = [cos(deg2rad(IMGdataset.Ki67.zRotation)), -sin(deg2rad(IMGdataset.Ki67.zRotation));...
     sin(deg2rad(IMGdataset.Ki67.zRotation)), cos(deg2rad(IMGdataset.Ki67.zRotation))];
 
IMGdataset.Ki67.ROIrotINV = [cos(deg2rad(-IMGdataset.Ki67.zRotation)), -sin(deg2rad(-IMGdataset.Ki67.zRotation));...
     sin(deg2rad(-IMGdataset.Ki67.zRotation)), cos(deg2rad(-IMGdataset.Ki67.zRotation))]; 
 
% register CD34 and CD3  
picMov = imread(IMGdataset.CD03.thumbfile);
[DataReg] = regThumbnails(picFix, picMov, 'Off');

IMGdataset.CD03.xMove = DataReg.xMove;
    IMGdataset.CD03.ROIxMove = - DataReg.xMove;
IMGdataset.CD03.yMove = DataReg.yMove;
    IMGdataset.CD03.ROIyMove = - DataReg.yMove;
IMGdataset.CD03.zRotation = DataReg.zRotation;

IMGdataset.CD03.thumfileReg = [cnst.regOutputFolder, IMGdataset.ID,'-CD03-ThumbReg.tif' ];
imwrite(DataReg.Thumb, IMGdataset.CD03.thumfileReg);
IMGdataset.CD03.ROIrot = [cos(deg2rad(IMGdataset.CD03.zRotation)), -sin(deg2rad(IMGdataset.CD03.zRotation));...
     sin(deg2rad(IMGdataset.CD03.zRotation)), cos(deg2rad(IMGdataset.CD03.zRotation))];

IMGdataset.CD03.ROIrotINV = [cos(deg2rad(-IMGdataset.CD03.zRotation)), -sin(deg2rad(-IMGdataset.CD03.zRotation));...
     sin(deg2rad(-IMGdataset.CD03.zRotation)), cos(deg2rad(-IMGdataset.CD03.zRotation))];
  
end % end of CalcRoiOffset