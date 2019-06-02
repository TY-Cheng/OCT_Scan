% 
%     {{Caserel}}
%     Copyright (C) {{2013}}  {{Pangyu Teng}}
% 
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 2 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License along
%     with this program; if not, write to the Free Software Foundation, Inc.,
%     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
%


function paths = getHyperReflectiveLayers(inputImg, constants, membraneIntensity)
% isos: membrane upper boundary

if nargin < 1
    display('requires at least 1 input (findHyperReflectiveZones.m)');
    return;
end

if nargin == 1
    %initiate parameters
    constants.shrinkScale = 0.3;
    constants.offsets = -3:3;
end

isPlot = 0;

%shrink the image.
szImg = size(inputImg);
procImg = imresize(inputImg,constants.shrinkScale,'bicubic');

%create adjacency matrices
[adjMatrixW, adjMatrixMW, adjMX, adjMY, adjMW, adjMmW, newImg] = getAdjacencyMatrix(procImg);

%create roi for getting shortestest path based on gradient-Y image.
[gx, gy] = gradient(newImg);
gyMinus = gy*-1; 
szImgNew = size(newImg);
roiImgmW = zeros(szImgNew);
roiImgmW(gy >= 0.8*mean(gy(:))) = 1 ;
roiImgW = zeros(szImgNew);
roiImgW(gyMinus >= mean(gyMinus(:))) = 1 ;

% find at least 2 layers
path{1} = 1;
count = 1;
while ~isempty(path) && count <= 2

    if count == 1
        % dark-to-light adjacency matrix
        roiImg = roiImgmW;
        adjM = adjMmW;
    else
        if membraneIntensity == 0
            % light-to-dark adjacency matrix
            roiImg = roiImgW;
            adjM = adjMW;
        end
    end
        
    %add columns of one at both ends of images
    roiImg(:,1)=1; roiImg(:,end)=1;
    
    includeX = ismember(adjMX, find(roiImg(:) == 1));
    includeY = ismember(adjMY, find(roiImg(:) == 1));
    % include only region of interst in the adjacency matrix
    keepInd = includeX & includeY;
    
    % compile adjacency matrix
    adjMatrix = sparse(adjMX(keepInd),adjMY(keepInd),adjM(keepInd),numel(newImg(:)),numel(newImg(:)));
    
    % get layer going from dark to light        
    [ dist,path{1} ] = graphshortestpath( adjMatrix, 1, numel(newImg(:)));
    
    if ~isempty(path{1})
                        
        % get rid of first few points and last few points
        [pathX,pathY] = ind2sub(szImgNew,path{1});        

        pathX = pathX(gradient(pathY)~=0);
        pathY = pathY(gradient(pathY)~=0);
        
        %block the obtained path and abit around it
        pathXArr = repmat(pathX,numel(constants.offsets),1);
        pathYArr = repmat(pathY,numel(constants.offsets),1);
        for i = 1:numel(constants.offsets)
            pathXArr(i,:) = pathXArr(i,:)+constants.offsets(i);
        end
        
        pathXArr = pathXArr(pathXArr > 0 & pathXArr <= szImgNew(1));
        pathYArr = pathYArr(pathXArr > 0 & pathXArr <= szImgNew(1));
        
        pathArr = sub2ind(szImgNew,pathXArr,pathYArr);
        if membraneIntensity == 0
            roiImgW(pathArr) = 0;
        else
            roiImg(pathArr) = 0;
        end
        paths(count).pathX = pathX;
        paths(count).pathY = pathY;

        if isPlot
            subplot(1,3,1);
            imagesc(procImg);
            subplot(1,3,2);
            imagesc(gy);        
            subplot(1,3,3);
            imagesc(roiImg);
            drawnow;
%             pause;
        end
        
    end % of ~empty
    count = count + 1;
end


if ~exist('paths','var')
    paths = {};
%     keyboard;
    return;
end % if exist

%format paths back to original size
for i = 1:numel(paths)    
    [paths(i).path, paths(i).pathY, paths(i).pathX] = resizePath(szImg, szImgNew, constants, paths(i).pathY, paths(i).pathX);    
    paths(i).pathXmean = nanmean(paths(i).pathX);
    paths(i).name = [];
    
end

%name each path (numel(paths) should equal to 2)
if numel(paths) ~= 2
    paths = {};
    display('error');
    return;
end

%based on the mean location detemine the layer type.
if paths(1).pathXmean < paths(2).pathXmean
    paths(1).name = 'ilm';
    paths(2).name = 'isos';
else
    paths(1).name = 'isos';    
    paths(2).name = 'ilm';    
end


if isPlot;
    figure;
    imagesc(inputImg);
    axis image; colormap('gray');
    hold on;
    for i = 1:numel(paths)
        cola = rand(1,3);
        plot(paths(i).pathY,paths(i).pathX,'r-','linewidth',1.5);
        text(paths(i).pathY(end),paths(i).pathX(end)-15,paths(i).name,'color',rand(1,3));
        drawnow;
    end
    hold off;
end

