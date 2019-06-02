excelCompiled = {};
excelLumped = {};

% intitialize a vector location of the layers
if params.membraneIntensity == 0
    layersToPlot  = {'ilm' 'isos' 'rpe'};
    % select the upper and lower boundaries of the interested layer
    layersToAnalyze = {'ilm' 'isos'};
end

% format the paths for analysis (get rid of unneeded 
imageLayer = formatPathsForAnalysis(imageLayer);

%% iterate through 'imageLayer(i).retinalLayers(j)'
% and save location to the corresponding vector 'layerCompile(storeInd)'

% quantify layer thickness
excel = {};

    if params.membraneIntensity == 1
        if numel(imageLayer.layers) == 3
            layersToPlot  = {'inlopl' 'ilm' 'isos'};
            % select the upper and lower boundaries of the interested layer
            layersToAnalyze = {'ilm' 'isos'};
        else
            layersToPlot  = {'ilm' 'isos'};
            % select the upper and lower boundaries of the interested layer
            layersToAnalyze = {'ilm' 'isos'};
        end
    end
    
    for k = 1:numel(layersToPlot)
        layerCompile(k).name = layersToPlot{k};
        layerCompile(k).x = [];
    end
    
    for j = 1:numel(imageLayer.layers)

        %find location in layerCompile to save the new pathX
        storeInd = find( strcmpi(imageLayer.layers(j).name,layersToPlot) ==1);

        if ~isempty(storeInd)
            layerCompile(storeInd).x = [ layerCompile(storeInd).x imageLayer.layers(j).pathXAnalysis];
        end

    end % of for j = 1:numel(imageLayer(i).retinalLayers),

    for k = 2:numel(layersToAnalyze)
        firstLayerInd = find(strcmpi(layersToAnalyze{k-1},layersToPlot)==1);
        secondLayerInd = find(strcmpi(layersToAnalyze{k},layersToPlot)==1);
        
        secondLayerx = layerCompile(secondLayerInd).x;
        firstLayerx = layerCompile(firstLayerInd).x;
        Z = secondLayerx - firstLayerx;
        Zmean = nanmean(Z);
        Ra = nanmean(abs(Z-Zmean));
        Ra1 = nanmean(abs(Z-Zmean)/Zmean);
        
        excel = [excel; { imageLayer.filename, strcat( [ layersToAnalyze{k-1} ' - ' layersToAnalyze{k}] ),...
            Zmean,...
            nanstd(Z),...
            Ra, Ra1}];
    end



% write in txt file
fprintf(fileID, '%s,%s,%5.5f,%5.5f,%5.5f,%5.5f\r\n', imageLayer.filename, ...
    strcat( [ layersToAnalyze{k-1} ' - ' layersToAnalyze{k}] ), ...
    Zmean, ...
    nanstd(Z), ...
    Ra, ...
    Ra1);
% T = cell2table(excel, 'VariableNames', {'image' 'layer' 'mean' 'sd' 'roughness_Ra' 'roughness_Ra1'});
% % writetable(T, [folderPath ,'./',folder,'.txt']);

% print out thickness
% T


