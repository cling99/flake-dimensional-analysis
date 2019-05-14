function afmImgInfo = getAFMFlakes(file, contrast, edgeMethod, fudgeFactor,...
                                    filterRadius, dilateRadius)

% Adds IBW files to path, 
% supresses warnings from APP version which stores these file in
% different locations
warning('off', 'MATLAB:mpath:nameNonexistentOrNotADirectory');
addpath('.\IBW', '-frozen');
warning('on', 'MATLAB:mpath:nameNonexistentOrNotADirectory');

ibwRaw = IBWread(file);

% For calculating flake area and dimensions
pixelX = ibwRaw.dx(1);
pixelY = ibwRaw.dx(2);
pixelArea = pixelX * pixelY;

amplitudeData = ibwRaw.y(:,:,1);
[BWfinal, Segout, overlay] = imgProcessFlakes();
%figure, imshow(BWfinal), title('segmented image');
%figure, imshow(Segout), title('outlined original image');

% For Labeling Flakes, Converts to RGB so we can write directly on img
h = figure;
set(h,'visible','off'), imshow(overlay), title('outlined original image (labeled)');
% Increase size
% set(h,'Position',[10 10 1000 1000])

% Finds and label connected components (starts at 1)
[connRegLabel, connRegNum] = bwlabel(BWfinal);

%% Iterates through different detected flakes on the image
% Since flakes might be rejected, flakeInfo array unallocated
flakeInfo = [];
flakeNum = 0;
for ii = 1:connRegNum
    flakeLogicArray = (connRegLabel == ii);
    
    % If "Flake" is greater than a single spot
    if nnz(flakeLogicArray) > 9
        % Finds flake area by summing across the number of pixels considered to
        % be part of the flake
        flakeArea = pixelArea * sum(sum(flakeLogicArray));
        
        % Finds the max height and mean height
        correctedAmpData = correctAmplitude2(amplitudeData, flakeLogicArray);
        
        maxHeight = max(max(correctedAmpData(flakeLogicArray)));
        % Converts array to vector first
        medianHeight = median(reshape(correctedAmpData(flakeLogicArray),...
                                      [], 1));
        meanHeight = mean(mean(correctedAmpData(flakeLogicArray)));
        
        % Finds Flake Volume
        flakeVol = sum(sum(1e9 * correctedAmpData(flakeLogicArray)...
                    .* (pixelArea * (1e9)^2)));
        
        % Approximates a length using region props
        [flakePosRow, flakePosCol] = find(flakeLogicArray);
        flakeLength = findMaxLength([flakePosRow, flakePosCol]);
        
        % Marks Flake on image
        flakeNum = flakeNum + 1;
        flakeStrNum = num2str(flakeNum);
        midEl = floor(length(flakePosCol)/2);
        
        % Inserts labels to labeled image
        text(flakePosCol(midEl), flakePosRow(midEl), flakeStrNum,...
                        'color', 'r',...
                        'fontsize', floor(min(size(BWfinal)) / 30));
                        % Alters font size based on image size
        
        flakeInfo(end+1,:) = [flakeArea, maxHeight, meanHeight,...
                                flakeLength, flakeVol, medianHeight]';
    end
end

% Take a screenshot of the image, rotate to match AFM
camroll(90); labeledImgNum = getframe;
close(1);

flakeInfo(:,[2:4 6]) = flakeInfo(:,[2:4, 6]) .* 1e9; % Scales to nm
flakeInfo(:,1) = flakeInfo(:,1) .* (1e9)^2; % Scales to nm^2

afmImgInfo.area = flakeInfo(:,1);
afmImgInfo.length = flakeInfo(:,4);
afmImgInfo.meanThick = flakeInfo(:,2);
afmImgInfo.medianThick = flakeInfo(:,6);
afmImgInfo.maxThick = flakeInfo(:,3);
afmImgInfo.volume = flakeInfo(:,5);
afmImgInfo.rawAFMData = amplitudeData;
% Rotate to match with Gwyddion
afmImgInfo.labeledFlakes = imrotate(Segout, 90); 
afmImgInfo.labeledFlakesNum = labeledImgNum.cdata;

afmImgInfo.flakeInfo = flakeInfo;
                   
%% Uses matlab's image processing toolbox to get flake locations
    function [BWfinal, Segout, overlay] = imgProcessFlakes()
        %figure, surf(amplitudeData), title('Raw Amplitude Data');
        
        img = mat2gray(amplitudeData);
        
        % Increase contrast if option selected
        if contrast
            img = imadjust(img); % Increase contrast
        end
        %figure, imshow(img), title('Unfiltered Image');
        img = wiener2(img, [filterRadius filterRadius]);
        %figure, imshow(img), title('Adaptive filtered Image');
        
        % Code from matlab example for detecting cell
        % https://www.mathworks.com/help/images/detecting-a-cell-using-image-segmentation.html
        
        [~, threshold] = edge(img, edgeMethod);
        BWs = edge(img,edgeMethod, threshold * fudgeFactor);
        %figure, imshow(BWs), title('binary gradient mask');
        
        % Dilating the image
        se90 = strel('line', dilateRadius, 90);
        se0 = strel('line', dilateRadius, 0);
        
        BWsdil = imdilate(BWs, [se90 se0]);
        %figure, imshow(BWsdil), title('dilated gradient mask');
        
        % Fill Interior Gaps
        BWdfill = imfill(BWsdil, 'holes');
        %figure, imshow(BWdfill);
        %title('binary image with filled holes');
        
        % Removing objects on border
        BWnobord = imclearborder(BWdfill, 4);
        %figure, imshow(BWnobord), title('cleared border image');
        
        % Smoothen the objects
        seD = strel('diamond',1);
        BWfinal = imerode(BWnobord,seD);
        BWfinal = imerode(BWfinal,seD);
        %figure, imshow(BWfinal), title('segmented image');
        
        % Applies various masks for output
        BWoutline = bwperim(BWfinal);
        SegoutR = img;
        SegoutG = img;
        SegoutB = img;
        
        % Converts Segout from white to Red
        SegoutR(BWoutline) = 255;
        SegoutG(BWoutline) = 0;
        SegoutB(BWoutline) = 0;
        Segout = cat(3, SegoutR, SegoutG, SegoutB);
        
        %figure, imshow(Segout), title('outlined original image');
        
        % Label Overlay
        overlay = labeloverlay(img, BWfinal);
    end

%% Finds the max length of a flake given a set of points
% INPUTS    POS: All positions of the flake, as a list of position vectors
% OUTPUTS   LENGTH: The maximum length of the flake in mm
function maxLength = findMaxLength(pos)
% If only 1 pixel (1 row), NaN
numPixels = size(pos,1);

maxLength = -inf;
% Iterates through each position combos
for ii = 1:numPixels
    for jj = ii:numPixels
        distance = sqrt((pixelX * (pos(ii,1) - pos(jj,1))).^2 +...
            (pixelY * (pos(ii,2) - pos(jj,2))).^2 );
        if distance > maxLength
            maxLength = distance;
        end
    end
end
end
end

%% Returns amplitude array scale by lowest element in filter
%   INPUTS  A: 2D matrix to create linear base
%           filter: A logic matrix returned by BWLabel
%   OUTPUTS correctA: A with corrected height around filter
%           a, b, c: coefficient for vector multiplication
function correctedMat = correctAmplitude2(A, filter)
minVal = min(min(A(filter)));
correctedMat = A - minVal;
end

%% Converts an binary BW image to an RGB image
% INPUTS    A binary BW image
% OUTPUTS   An rgb image
function rgb = bw2rgb(bw)
    grayImage = 255 * uint8(bw);
    rgb = cat(3, grayImage, grayImage, grayImage);
end
