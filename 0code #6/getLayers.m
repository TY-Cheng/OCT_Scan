function [layers, params] = getLayers(img,params)
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

%initialize constants
if nargin < 2        
    
    % resize the image if 1st value set to 'true',
    % with the second value to be the scale.
    params.isResize = [true 0.5];     
        
    % constants used for defining the region for segmentation of individual layer
    params.roughILMandISOS.shrinkScale = 0.2;
    params.roughILMandISOS.offsets = -20:20;    
    params.ilm_0 = 4;
    params.ilm_1 = 4;
    params.isos_0 = 4;
    params.isos_1 = 4;
    params.rpe_0 = 0.1;
    params.rpe_1 = 0.01;
    params.inlopl_0 = 0.4; %   0.4;% 0.1
    params.inlopl_1 = 0.5; %   0.5;% 0.3
    params.nflgcl_0 = 0.01;%  0.01;  0.05
    params.nflgcl_1 = 0.1; %   0.1;  0.3
    params.iplinl_0 = 0.6;
    params.iplinl_1 = 0.2;
    params.oplonl_0 = 0.05;%4;
    params.oplonl_1 = 0.5;%4;    
        
    % parameters for ploting
    params.txtOffset = -7;
    colorarr=colormap('jet'); 
    params.colorarr=colorarr(64:-8:1,:);
    
    % a constant (not used in this function, used in 'octSegmentationGUI.m'.)
    params.smallIncre = 2;    
    
end

%clear up matlab's mind
clear retinalLayers

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
% layerSegmentationOrder = {'roughILMandISOS' 'ilm' 'isos' 'rpe' 'inlopl' 'nflgcl' 'iplinl' 'oplonl'};
layerSegmentationOrder = {'roughILMandISOS' 'ilm' 'isos' 'rpe'};
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

%format paths back to original size
szImgNew = size(imgNew);
constants.shrinkScale = params.isResize(2);
for i = 1:numel(layers)    
    [layers(i).path, layers(i).pathY, layers(i).pathX] = resizePath(szImg, szImgNew, constants, layers(i).pathY, layers(i).pathX);    
    layers(i).pathXmean = nanmean(layers(i).pathX);
    
end

% plot oct image and the obtained retinal layers.
isPlot = 1;
if isPlot

    imagesc(img_full);
    axis image; colormap('gray'); hold on; drawnow;

    layersToPlot = {'ilm' 'isos' 'rpe'};
    hOffset =       [40    0      40    0        0        40       -40      -40]; % for displaying text
    for k = 1:numel(layersToPlot)

        matchedLayers = strcmpi(layersToPlot{k},{layers(:).name});
        layerToPlotInd = find(matchedLayers == 1);

        if ~isempty(layers(layerToPlotInd).pathX)            
            colora = params.colorarr(k,:);
            plot(layers(layerToPlotInd).pathY,layers(layerToPlotInd).pathX-1,'color',colora,'linewidth',1.5);
            plotInd = round(numel(layers(layerToPlotInd).pathX)/2);            
            text(layers(layerToPlotInd).pathY(plotInd)+hOffset(k),layers(layerToPlotInd).pathX(plotInd)+params.txtOffset,layers(layerToPlotInd).name,'color',colora,'linewidth',2);            
            drawnow;
        end % of if ~isempty            

    end % of k
    hold off;

end % of isPlot        
