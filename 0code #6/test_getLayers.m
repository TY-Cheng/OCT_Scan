close all;clear all;clc;
%% Section 0, user inputs

% IMPORTANT: select 0 for membrane darker than biomass, or 1 otherwise
params.membraneIntensity = 1;

% choose the folder containing tif images
folder = '1';

% select 1 to plot segmented images, or 0 not to plot
params.isPlot = 0;

% select 1 to crop images, or 0 for no cropping
cropimages = 0;

%% Section 1, loads the path of the image.

path = ['./',folder,'/'];
[filename, folderPath , filterindex] = uigetfile([path '*.tif' ],'Pick some images','MultiSelect', 'on');
if ~iscell(filename)
    filename = {filename};
end
for i = 1:numel(filename)
    imagePath{i} = [folderPath ,filename{i}];
end

yrange = []; xrange = [];
if cropimages == 1
    figure;
    title('pick a region of interest to segment for the selected images');
    [trsh, rect] = imcrop(imread(imagePath{1}));
    xrange = round(rect(1)):round(rect(1)+rect(3));
    yrange = round(rect(2)):round(rect(2)+rect(4));
end

%% Section 2, automatically segments the 2 layers based on graph theory.

fileID = fopen([folderPath ,'./',folder,'.txt'],'a');
fprintf(fileID, '%s,%s,%s,%s,%s,%s\r\n', 'image','layer','mean','sd','roughness_Ra','roughness_Ra1');

for i = 1:numel(imagePath)
    tic
    try
        clearvars -except fileID filename  folderPath i imagePath  path
        params.membraneIntensity = 1;   folder = '1';   params.isPlot = 0;
        yrange = []; xrange = [];
        %%
        % read in the image.
        img = imread(imagePath{i});
        
        % error checking, get one channel from image.
        if size(img,3) > 1
            img = img(:,:,1);
            fprintf('warning: this is probably not an oct image');
        end
        
        % make image type as double.
        img = double(img);
        
        % get size of image.
        szImg = size(img);
        
        %segment whole image if yrange/xrange is not specified.
        if isempty(yrange) && isempty(xrange)
            yrange = 1:szImg(1);
            xrange = 1:szImg(2);
        end
        img = img(yrange,xrange);
        
        %%%------------------- our denoising method ------------------------%%%
        fprintf('denoising image %d of %d\n',i,numel(imagePath));
        k = 1.8; % distribution coefficient estimation: k ~= mean/std
        c1 = (1-1/(2*k*k))^(1/4);
        c2 = 1-(1-1/(2*k*k))^(1/2);
        
        % parameter selection
        par.gamma = 10;
        par.tau = 0.9/(1+par.gamma*(2^2));
        par.lambda = 0.2; % regularizer weight
        par.theta = 1;
        par.c1 = c1;
        par.c2 = c2;
        par.maxIter = 15;
        
        % modify input to correct range
        D_input = ( img-min(img(:)) )/( max(img(:))-min(img(:)) );
        [ img_clean ] = ladexp_huberTV( D_input, 1, par );
        %     img_clean = img_clean.^(1.25);
        %     img_clean = 255*(img_clean - min(img_clean(:)))/( max(img_clean(:)) - min(img_clean(:)) );
        
        %%%----------------- contrast and brightness -----------------------%%%
        img_clean = imadjust(img_clean,[min(img_clean(:)); max(img_clean(:))],[0; 1], 2.5);
        %     figure(2);clf; imshow([D_input;img_clean],[]);
        
        %     % parameter for smothing the images.
        %     params.filter0Params = [5 5 1];
        %     %smooth image with specified kernels
        %     %for denosing
        %     img_clean = imfilter(img_clean,fspecial('gaussian',params.filter0Params(1:2),params.filter0Params(3)),'replicate');
        
        
        fprintf('segmenting image %d of %d\n',i,numel(imagePath));
        % get layers for different cases
        if params.membraneIntensity == 0
            [layers, params] = getLayers_darkMembrane(img_clean, params);
        else
            [layers, params] = getLayers_brightMembrane(img_clean, params);
        end
        % save range of image.
        params.yrange = yrange;
        params.xrange = xrange;
        
        % save data to struct.
        imageLayer.filename = filename{i};
        imageLayer.imagePath = imagePath{i};
        imageLayer.layers = layers;
        imageLayer.params = params;
        
        % save segmentation
        filepath = [imageLayer.imagePath '_octSegmentation.mat'];
        save(filepath, 'imageLayer');
        fprintf('segmentation saved to %s\n',filepath);
        
        
        %%  Section 4, calculate and print out retinal thickness (in pixels)
        
        calculateRetinalThickness
    catch
        fprintf(['\nCannot read:\n', imagePath{i}, '\n\n'])
    end
    toc
end

fclose(fileID);


%%   Section 3, using a GUI, iterate through the segmentation results,
%              and maually or semi-automatically correct the segmented
%              layers.
%
% close all;
%
% filename = [imagePath{1} '_octSegmentation.mat'];
%
% isReviewSegmentation = 1;
% if isReviewSegmentation,
%     [h,guiParam] = octSegmentationGUI(filename);
%
%     if guiParam.proceed
%         delete(guiParam.figureOCT);
%         delete(h);
%     else
%         return;
%     end
% end


