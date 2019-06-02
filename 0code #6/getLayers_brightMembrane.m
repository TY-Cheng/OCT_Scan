function [layers, params] = getLayers_brightMembrane(img,params)
%%
%
% [layers params] = getLayers(img,params)
% identifies the boundaries between layeres given an optical 
% coherence tomography image, 'img'.
% 
% The method for identification of these layer boundaries is
% based on graph theory.
% 
%%

if nargin < 1
    display('requires 1 input');
    return;
end

%initialize default parameters      
    
    % resize the image if 1st value set to 'true',
    % with the second value to be the scale.
    params0.isResize = [true 0.5];     
        
    % constants used for defining the region for segmentation of individual layer
    params0.roughILMandISOS.shrinkScale = 0.5;
    params0.roughILMandISOS.offsets = -1:1;  
    params0.membraneIntensity = 1;
    params0.ilm_0 = 8;
    params0.ilm_1 = 100;
    params0.isos_0 = 2;
    params0.isos_1 = 50;
    params0.rpe_0 = 0.1;
    params0.rpe_1 = 0.01;
    params0.inlopl_0 = 0.1; %   0.4;% 0.1
    params0.inlopl_1 = 0.08; %   0.5;% 0.3
    params0.nflgcl_0 = 0.01;%  0.01;  0.05
    params0.nflgcl_1 = 0.1; %   0.1;  0.3
    params0.iplinl_0 = 0.6;
    params0.iplinl_1 = 0.2;
    params0.oplonl_0 = 0.05;%4;
    params0.oplonl_1 = 0.5;%4;    
        
    % parameters for ploting
    params0.txtOffset = -7;
    colorarr=colormap('jet'); 
    params0.colorarr=colorarr(64:-8:1,:);
    
    % a constant (not used in this function, used in 'octSegmentationGUI.m'.)
    params0.smallIncre = 2;    
    params0.isPlot = 0;

    if exist('params','var')
        %merge parameters
        params = mergeParam(params0, params);
    else
        params = params0;
    end

%clear up matlab's mind
clear layers paths

%get image size
szImg = size(img);
img_full = img;

%resize image.
if params.isResize(1)
    img = imresize(img,params.isResize(2),'bicubic');
end

% create adjacency matrices and its elements base on the image.
[params.adjMatrixW, params.adjMatrixMW, params.adjMA, params.adjMB, params.adjMW, params.adjMmW, imgNew] = getAdjacencyMatrix(img);

% obtain rough segmentation of the ilm and isos, then find the 
% layers in the order of 'layerSegmentationOrder'
%%vvvvvvvvvvvvvvvDO  NOT  CHANGE BELOW LINE (ORDER OF LAYERS SHALL NOT BE CHANGED)vvvvvvvvvvvvvv%%
% isos: membrane upper boundary
% layerSegmentationOrder = {'roughILMandISOS' 'ilm' 'isos' 'inlopl'};
layerSegmentationOrder = {'roughILMandISOS' 'isos' 'ilm' 'inlopl'};
%%^^^^^^^^^^^^^^^DO  NOT  CHANGE ABOVE LINE (ORDER OF LAYERS SHOULD NOT BE CHANGED)^^^^^^^^^^^^^%%
% segment retinal layers
layers = [];
for layerInd = 1:numel(layerSegmentationOrder)        
    [layers, ~] = getLayersCore(layerSegmentationOrder{layerInd},imgNew,params,layers);
end

%delete elements of the adjacency matrices prior function exit to save memory
toBeDeleted = {'adjMatrixWSmo' 'adjMatrixMWSmo' 'adjMWSmo' 'adjMmWSmo'  'adjMW' 'adjMmW' 'adjMatrixW' 'adjMatrixMW' 'adjMA' 'adjMB'};
for delInd = 1:numel(toBeDeleted)
    params.(toBeDeleted{delInd}) = [];
end

% plot oct image and the obtained layers.
if params.isPlot
    paths = [];
    szImgNew = size(imgNew);
    constants.shrinkScale = params.isResize(2);
    %format paths back to original size
    for i = 1:numel(layers)    
        [paths(i).path, paths(i).pathY, paths(i).pathX] = resizePath(szImg, szImgNew, constants, layers(i).pathY, layers(i).pathX);    
        paths(i).pathXmean = nanmean(layers(i).pathX);
        paths(i).name = layers(i).name;
    end
    
    figure(1);clf;imagesc(img_full);
    axis image; colormap('gray'); hold on; drawnow;
    
    if layers(3).dist > (0.6*2)^1*szImgNew(2) % 0.5*(2*numel(layers(3).pathY))
        layers(3) = [];
        layersToPlot = {'ilm' 'isos'};
    else
        layersToPlot = {'ilm' 'isos' 'inlopl'};
    end
    hOffset =       [40    0      40    0        0        40       -40      -40]; % for displaying text
    for k = 1:numel(layersToPlot)

        matchedLayers = strcmpi(layersToPlot{k},{layers(:).name});
        layerToPlotInd = find(matchedLayers == 1);

        if ~isempty(layers(layerToPlotInd).pathX)            
            colora = params.colorarr(k,:);
            plot(paths(layerToPlotInd).pathY,paths(layerToPlotInd).pathX-1,'--','color',colora,'linewidth',1.5);
            plotInd = round(numel(paths(layerToPlotInd).pathX)/2);            
            text(paths(layerToPlotInd).pathY(plotInd)+hOffset(k),paths(layerToPlotInd).pathX(plotInd)+params.txtOffset,paths(layerToPlotInd).name,'color',colora,'linewidth',2);            
            drawnow;
        end % of if ~isempty            

    end % of k
    hold off;

end % of isPlot     
