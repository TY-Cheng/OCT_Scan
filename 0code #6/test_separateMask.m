close all;clear all;clc;
%% Section 0, user inputs

    % choose the folder containing tif mask images
    folder = 'mask_images'; 
    
    % select 1 to plot segmented images, or 0 not to plot
    params.isPlot = 1;

    % select 1 to crop images, or 0 for no cropping
    cropimages = 0;
    
    SaveFolderName = folder;
    mkdir('results',SaveFolderName);

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
        [trsh rect] = imcrop(imread(imagePath{1}));
        xrange = round(rect(1)):round(rect(1)+rect(3));
        yrange = round(rect(2)):round(rect(2)+rect(4));
    end
    
%% Section 2, automatically segments the 2 layers based on graph theory.

for i = 1:numel(imagePath)
    
    % read in the image.
    InfoImage=imfinfo(imagePath{i});
    mImage=InfoImage(1).Width;
    nImage=InfoImage(1).Height;
    NumberImages=length(InfoImage);
 
    img=zeros(nImage,mImage,NumberImages,'uint16');
    for j=1:NumberImages
       img(:,:,j)=imread(imagePath{i},'Index',j);
    end
    
    % make image type as double.
    img = double(img);
%     img = img(30:end,:,:);
    
    % get size of image.
    szImg = size(img);
    maskBiofilm=zeros(szImg);
    
    %segment whole image if yrange/xrange is not specified.
    if isempty(yrange) && isempty(xrange)
        yrange = 1:szImg(1);
        xrange = 1:szImg(2);
    end    
    img = img(yrange,xrange,:);   
    % save range of image.
    params.yrange = yrange;
    params.xrange = xrange;
    
    % modify input to correct range
    maskIm = ( img-min(img(:)) )/( max(img(:))-min(img(:)) );
    
    %%%-----------------------------------------------------------------%%%
    
    for j=1:NumberImages
        fprintf('segmenting image %d of %d\n',j,numel(imagePath));
        % extract biofilm mask part
        [maskBiofilm(:,:,j), params] = getMask_biofilm(maskIm(:,:,j), params);
    end
    
    saveastiff(uint8(255*maskBiofilm), fullfile('results',SaveFolderName,[filename{i}(1:end-4),'_biofilm.tif']));
    
end
