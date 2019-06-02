function [rPaths, img] = getBiofilmCore(layerName,img,params,rPaths)

% this is a founction used within getRetinalLayers.m

if nargin < 3
    display('3 inputs required, getLayers.m');
    return;   
end

szImg = size(img);

switch layerName
    
    case {'roughILMandISOS'}
        
        imgOld = img(:,2:end-1);
        pathsTemp = getHyperReflectiveLayers_Biofilm(imgOld,params.roughILMandISOS,params.membraneIntensity);                 
                        
        %save to structure 
        clear rPaths
        rPaths = pathsTemp;
        
        return;        

    case {'nflgcl' 'inlopl' 'ilm' 'isos' 'oplonl' 'iplinl' 'rpe'}
        
        adjMA = params.adjMA;
        adjMB = params.adjMB;
        adjMW = params.adjMW;
        adjMmW = params.adjMmW;        
        
    case {'IF_YOU_WANT_A_SMOOTHER_EDGE_PUT_THE_LAYER_NAMES_HERE'}

        adjMA = params.adjMA;
        adjMB = params.adjMB;
        adjMW = params.adjMWSmo;
        adjMmW = params.adjMmWSmo;
        
end


% initialize region of interest
szImg = size(img);
roiImg = zeros(szImg);

% avoid the top part of image
roiImg(1:20,:) = 0;

% select region of interest based on layers priorly segmented.
for k = 2:szImg(2)-1
    
    switch layerName
        
        case {'nflgcl'}
            
            % define a region (from 'startInd' to 'endInd') between 'ilm'
            % and 'inlopl'.
            indPathX = find(rPaths(strcmp('ilm',{rPaths.name})).pathY==k);
            startInd0 = rPaths(strcmp('ilm',{rPaths.name})).pathX(indPathX(1));            
            indPathX = find(rPaths(strcmp('inlopl',{rPaths.name})).pathY==k);
            endInd0 = rPaths(strcmp('inlopl',{rPaths.name})).pathX(indPathX(1));
                        
            startInd = startInd0 - ceil(params.nflgcl_0*(endInd0-startInd0));
            endInd = endInd0 - round(params.nflgcl_1*(endInd0-startInd0));
            
        case {'rpe'} % membrane lower boundary

            indPathX = find(rPaths(strcmp('isos',{rPaths.name})).pathY==k);
            
            % define a region (from 'startInd' to 'endInd') below 'isos'.
            startInd0 = rPaths(strcmp('isos',{rPaths.name})).pathX(indPathX(1));
            endInd0 = startInd0+max(round((rPaths(strcmp('isos',{rPaths.name})).pathXmean-rPaths(strcmp('ilm',{rPaths.name})).pathXmean)),30);

            startInd = startInd0+max(round(params.rpe_0*(endInd0-startInd0)),5);
            endInd = endInd0-round(params.rpe_1*(endInd0-startInd0));                   
            
        case {'inlopl'} % biomass rough upper boundary
            
            % define a region (from 'startInd' to 'endInd') above 'isos'.
%             indPathX = find(rPaths(strcmp('ilm',{rPaths.name})).pathY==k);
%             startInd0 = rPaths(strcmp('ilm',{rPaths.name})).pathX(indPathX(1));
%             indPathX = find(rPaths(strcmp('isos',{rPaths.name})).pathY==k);
%             endInd0 = rPaths(strcmp('isos',{rPaths.name})).pathX(indPathX(1));
            indPathX = find(rPaths(strcmp('ilm',{rPaths.name})).pathY==k);
            startInd0 = 10;
            endInd0 = rPaths(strcmp('ilm',{rPaths.name})).pathX(indPathX(1));
            
            startInd = startInd0+round(params.inlopl_0*(endInd0-startInd0));
            endInd = endInd0-round(params.inlopl_1*(endInd0-startInd0));            
            
        case {'ilm'} % biomass compact upper boundary
            
            % define a region (from 'startInd' to 'endInd') near 'ilm'.
            indPathX = find(rPaths(strcmp('ilm',{rPaths.name})).pathY==k);
            startInd = rPaths(strcmp('ilm',{rPaths.name})).pathX(indPathX(1)) - params.ilm_0; 
            endInd = rPaths(strcmp('ilm',{rPaths.name})).pathX(indPathX(1)) + params.ilm_1;             
            if startInd > rPaths(strcmp('ilm',{rPaths.name})).pathXmean - 6
                startInd = rPaths(strcmp('ilm',{rPaths.name})).pathXmean - 6;
            end
        case {'isos'} % membrane upper boundary  
            
            % define a region (from 'startInd' to 'endInd') near 'isos'.
            indPathX = find(rPaths(strcmp('isos',{rPaths.name})).pathY==k);            
            
            startInd = rPaths(strcmp('isos',{rPaths.name})).pathX(indPathX(1)) - params.isos_0; 
            endInd = rPaths(strcmp('isos',{rPaths.name})).pathX(indPathX(1)) + params.isos_1;             
            
        case {'iplinl'}

            % define a region (from 'startInd' to 'endInd') between
            % 'nflgcl' and 'inlopl'.
            indPathX = find(rPaths(strcmp('nflgcl',{rPaths.name})).pathY==k);
            startInd0 = rPaths(strcmp('nflgcl',{rPaths.name})).pathX(indPathX(1));
            indPathX = find(rPaths(strcmp('inlopl',{rPaths.name})).pathY==k);
            endInd0 = rPaths(strcmp('inlopl',{rPaths.name})).pathX(indPathX(1));
            
            startInd = startInd0 + round(params.iplinl_0*(endInd0-startInd0));
            endInd = endInd0 - round(params.iplinl_1*(endInd0-startInd0));
            
            
        case {'oplonl'}
                        
            % define a region (from 'startInd' to 'endInd') between
            % 'inlopl' and 'isos'.
            indPathX = find(rPaths(strcmp('inlopl',{rPaths.name})).pathY==k);
            startInd0 = rPaths(strcmp('inlopl',{rPaths.name})).pathX(indPathX(1));
            indPathX = find(rPaths(strcmp('isos',{rPaths.name})).pathY==k);
            endInd0 = rPaths(strcmp('isos',{rPaths.name})).pathX(indPathX(1));
                        
%           startInd = startInd0 + params.oplonl_0;
%           endInd = endInd0 - params.oplonl_1;
            startInd = startInd0 +round(params.oplonl_0*(endInd0-startInd0));
            endInd = endInd0 -round(params.oplonl_1*(endInd0-startInd0));
    end
    
    %error checking    
    if startInd > endInd
        startInd = endInd - 1;
    end            
    
    if startInd < 1
        startInd = 1;
    end
    
    if endInd > szImg(1)
        endInd = szImg(1);
    end
                    
    % set region of interest at column k from startInd to endInd
    roiImg(startInd:endInd,k) = 1;
    
end

%ensure the 1st and last column is part of the region of interest.
roiImg(:,1)=1;
roiImg(:,end)=1;            

% include only region of interst in the adjacency matrix
includeA = ismember(adjMA, find(roiImg(:) == 1));
includeB = ismember(adjMB, find(roiImg(:) == 1));
keepInd = includeA & includeB;

%     %alternative to ismember, 
%     roiImgOne = find(roiImg(:) == 1)';
%     includeA = sum(bsxfun(@eq,adjMA(:),roiImgOne),2);
%     includeB = sum(bsxfun(@eq,adjMB(:),roiImgOne),2);
%     keepInd = includeA & includeB;

%get the shortestpath
switch layerName
    %bright to dark
    case {'nflgcl' 'oplonl' 'iplinl' 'rpe'}
        adjMatrixW = sparse(adjMA(keepInd),adjMB(keepInd),adjMW(keepInd),numel(img(:)),numel(img(:)));    
        [ dist, path ] = graphshortestpath( adjMatrixW, 1, numel(img(:)) );
        % dist = nan(size(path));
        % for i = 1:numel(path)-1,dist(i)=adjMatrixW(path(i),path(i+1));end
    %dark to bright
    case {'inlopl' 'ilm'}
        adjMatrixMW = sparse(adjMA(keepInd),adjMB(keepInd),adjMmW(keepInd),numel(img(:)),numel(img(:)));    
        [ dist, path ] = graphshortestpath( adjMatrixMW, 1, numel(img(:)) );        
        % dist = nan(size(path));
        % for i = 1:numel(path)-1,dist(i)=adjMatrixMW(path(i),path(i+1));end
    %membrane dependent
    case {'isos'}
        if params.membraneIntensity == 0
            adjMatrixW = sparse(adjMA(keepInd),adjMB(keepInd),adjMW(keepInd),numel(img(:)),numel(img(:)));    
            [ dist, path ] = graphshortestpath( adjMatrixW, 1, numel(img(:)) );
        else
            adjMatrixMW = sparse(adjMA(keepInd),adjMB(keepInd),adjMmW(keepInd),numel(img(:)),numel(img(:)));    
            [ dist, path ] = graphshortestpath( adjMatrixMW, 1, numel(img(:)) );
        end
end

%convert path indices to subscript
[pathX, pathY] = ind2sub(szImg,path);

%if name layer existed, overwrite it, else add layer info to struct
matchedLayers = strcmpi(layerName,{rPaths(:).name});
layerToPlotInd = find(matchedLayers == 1);
if isempty(layerToPlotInd)    
    layerToPlotInd = numel(rPaths)+1;
    rPaths(layerToPlotInd).name = layerName;
end

% save data.
rPaths(layerToPlotInd).path = path;
% rPaths(layerToPlotInd).dist = dist;
rPaths(layerToPlotInd).pathX = pathX;
rPaths(layerToPlotInd).pathY = pathY;
rPaths(layerToPlotInd).dist = dist;
rPaths(layerToPlotInd).pathXmean = mean(rPaths(layerToPlotInd).pathX(gradient(rPaths(layerToPlotInd).pathY)~=0));

%create an additional smoother layer for rpe and isos
isSmooth = 0;
if isSmooth
    switch layerName
        case {'rpe' 'isos'}       

            %find lines where pathY is on the image
            rpePathInd = gradient(pathY) ~= 0;

            % fit line with cubic smoothing spline
            lambda = 1E-6; %small means really smooth
            pathXpoly = pathX;
            pathYpoly = pathY;
           
           if numel(pathX(rpePathInd))>=2 && numel(pathY(rpePathInd))==numel(pathX(rpePathInd))
               [pathXpoly(rpePathInd), ~] = csaps(pathY(rpePathInd),pathX(rpePathInd),...
                  lambda,pathY(rpePathInd));               
    %             [pathXpoly(rpePathInd), ~] = csaps_pt(pathY(rpePathInd),pathX(rpePathInd),...
    %                lambda,pathY(rpePathInd));               
           end
           
            %add layer info to struct
            %layerToPlotInd = numel(rPaths)+1;
            %rPaths(layerToPlotInd).name = 'rpeSmooth';
            %update rpw layer info to struct
            rPaths(layerToPlotInd).pathX = round(pathXpoly);
            rPaths(layerToPlotInd).pathY = round(pathYpoly);            
            rPaths(layerToPlotInd).path = sub2ind(szImg,rPaths(layerToPlotInd).pathX,rPaths(layerToPlotInd).pathY);
            rPaths(layerToPlotInd).pathXmean = mean(rPaths(layerToPlotInd).pathX(gradient(rPaths(layerToPlotInd).pathY)~=0)); 
           
        otherwise
    end
end