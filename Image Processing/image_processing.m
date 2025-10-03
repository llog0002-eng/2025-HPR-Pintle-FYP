clear all; close all; clc;

%% User settings
processGif = false;
processSprayAngle = true;
allFrames = false;
monoChannel = 1;  % 1=Red, 2=Green, 3=Blue
folder = '\\ad.monash.edu\home\User090\ppin0001\Desktop\Pintle Data\19_9_25\';
filePath = append(folder,'1_1_4_2_70bar.mraw');
backgroundPath = '\\ad.monash.edu\home\User090\ppin0001\Desktop\Pintle Data\19_9_25\1_1_3_1_60bar.mraw'; % DO NOT CHANGE
outputPath16 = append(folder, 'avgMonoSubFrame.tif');
outputPathBW = append(folder,'avgMonoSubFrameBW.tif');
gifPath = append(folder,'bitshift.gif');

frameDelay = 0.01;
visualFrames = 50;
chunkSize = 100;

sigmaSmooth = 1;
backgroundFrame = 1;
maxFrames = 1000;
bitshiftAmount = 2;
maxPixelValue = 4095;
pintleDiameter_mm = 25;
background_threshold = 0.009;

% Distances from the pintle tip for spray angle calculation
distances_mm = [10, 20];

% Flow detection settings
sigmaFactor = 0.008; % DO NOT CHANGE 
maxPixelValue = 4095;    % 12-bit camera max value

% Binarisation settings
windowSize = 101; %  DO NOT CHANGE 
kValue = 0.3; %  DO NOT CHANGE 

rowMid = 50;
epsVal = 1e-6;
minPixelSize = 60;

% Camera info
frameWidth = 1024; frameHeight = 640; numChannels = 3; bitDepth = 12;
thresholdDarkFraction = 0.08;

%% File info
fid = fopen(filePath,'r'); fseek(fid,0,'eof'); fileBytes = ftell(fid); fclose(fid);
bytesPerPixel = ceil(bitDepth/8);
bytesPerFrame = frameWidth * frameHeight * numChannels * bytesPerPixel;
numFrames = floor(fileBytes / bytesPerFrame);
numFrames = min(numFrames, maxFrames);
fprintf('Total frames used: %d\n', numFrames);

%% Load background
bgData = double(readmraw(backgroundPath, backgroundFrame));
bgMono = bgData(:,:,monoChannel);

if allFrames
    startFrame = 1;  % start from the first frame
    I = double(readmraw(filePath, startFrame));
    Ibright = uint16(I);
    Ibright = bitshift(Ibright, bitshiftAmount);
    Ibright(Ibright > maxPixelValue) = maxPixelValue;

    % Convert to grayscale
    Igray = uint8(255 * mat2gray(rgb2gray(Ibright)));

    h = imshow(Igray); % initial display
    title(sprintf('Frame %d', startFrame)); % initial frame title
    imwrite(Igray, gifPath, 'gif', 'LoopCount', Inf, 'DelayTime', frameDelay); % first frame

    for i = startFrame+1 : startFrame + visualFrames - 1
        I = double(readmraw(filePath, i));
        Ibright = uint16(I);
        Ibright = bitshift(Ibright, bitshiftAmount);
        Ibright(Ibright > maxPixelValue) = maxPixelValue;

        Igray = uint8(255 * mat2gray(rgb2gray(Ibright)));

        set(h, 'CData', Igray);
        title(sprintf('Frame %d', i)); % update title with current frame
        drawnow;
    end

    fprintf('Preview GIF saved to: %s\n', gifPath);
end

exactFlowFrame = detectFlowFrames(filePath, monoChannel, 1, maxFrames, sigmaFactor, maxPixelValue);

%% Display GIF (optional)
if processGif
    I = double(readmraw(filePath, exactFlowFrame));
    Ibright = uint16(I);
    Ibright = bitshift(Ibright, bitshiftAmount);
    Ibright(Ibright > maxPixelValue) = maxPixelValue;

    % Convert to grayscale for GIF
    Igray = uint8(255 * mat2gray(rgb2gray(Ibright)));
    figure;
    h = imshow(Igray); % initial display
    imwrite(Igray, gifPath, 'gif', 'LoopCount', Inf, 'DelayTime', frameDelay); % first frame

    for i = exactFlowFrame+1 : exactFlowFrame + visualFrames - 1
        I = double(readmraw(filePath, i));
        Ibright = uint16(I);
        Ibright = bitshift(Ibright, bitshiftAmount);
        Ibright(Ibright > maxPixelValue) = maxPixelValue;

        Igray = uint8(255 * mat2gray(rgb2gray(Ibright)));

        set(h, 'CData', Igray);
        drawnow;

        imwrite(Igray, gifPath, 'gif', 'WriteMode', 'append', 'DelayTime', frameDelay);
    end

    fprintf('GIF saved to: %s\n', gifPath);
end

%% Pintle detection
bwPost = bgMono < background_threshold*max(bgMono(:));
bwPost = imfill(bwPost,'holes');
bwPost = bwareaopen(bwPost,50);

rowData = bwPost(rowMid,:);
transitions = diff(rowData);
startIdx = find(transitions==1,1,'first');
endIdx = find(transitions==-1,1,'last');

if ~isempty(startIdx) && ~isempty(endIdx)
    centerX = round((endIdx-startIdx)/2 + startIdx);
    centerY = rowMid;
    pintleBox = [startIdx, rowMid-10, endIdx-startIdx, 20];
else
    centerX = size(bwPost,2)/2;
    centerY = rowMid;
    pintleBox = [centerX-1, centerY-1, 2, 2];
end
fprintf('Detected pintle center at (X,Y) = (%d,%d)\n', centerX, centerY);

colRange = startIdx:endIdx;
colData = any(bwPost(:,colRange),2);
yTop = find(colData,1,'first');
yBottom = find(colData,1,'last');
fprintf('Pintle vertical extent from Y = %d to Y = %d\n', yTop, yBottom);

pintleWidthPixels = endIdx - startIdx;
pixelSize_mm = pintleDiameter_mm / pintleWidthPixels;
fprintf('Pixel size: %.4f mm/pixel\n', pixelSize_mm);

%% Convert distances to pixels
distances_pix = round(distances_mm / pixelSize_mm);    
yDetect = yBottom + distances_pix;   

%% Frame averaging and background division
if processSprayAngle
    accumMono = zeros(frameHeight, frameWidth, 'double');

    for f = exactFlowFrame:chunkSize:numFrames
        endFrame = min(f + chunkSize - 1, numFrames);
        I_chunk = readmraw(filePath, [f endFrame]);
        monoChunk = double(I_chunk(:,:,monoChannel,:)) ./ (bgMono + epsVal);
        accumMono = accumMono + sum(monoChunk, 4);
    end

    avgMono = accumMono / numFrames;
    avgMonoCropped = avgMono(yBottom:end, :);

    BW_sauvola = sauvolaSingle(avgMonoCropped, [windowSize windowSize], kValue);
    BW_inv = ~BW_sauvola;
    BW_sauvola_inv = bwareaopen(BW_inv, 5,4);
    BW_sauvola = ~BW_sauvola_inv;

    figure('Name', 'Sauvola Mask Results', 'NumberTitle', 'off');
    subplot(1,2,1); imshow(BW_sauvola);  
    title(sprintf('Sauvola Mask w=%d k=%.2f', windowSize, kValue));
    subplot(1,2,2); imshow(avgMonoCropped); hold on;
    hOverlay = imshow(BW_sauvola);  
    set(hOverlay, 'AlphaData', 0.3);  
    hold off; title('Sauvola Overlay on Cropped Image');

    imwrite(BW_sauvola, 'avgMonoSauvolaMask_cropped.jpg');
    fprintf('Cropped Sauvola mask saved as avgMonoSauvolaMask_cropped.jpg\n');

    %% Edge detection
    edgesCanny = edge(BW_sauvola, 'Canny', [], 3);

    %% Clean small noise and connect edges
    edgesCannyClean = bwareaopen(edgesCanny, minPixelSize);
    se = strel('line', 5, 0);                       
    edgesCannyClean = imclose(edgesCannyClean, se);

    %% Export edge image
    edgeCannyPath = strrep(outputPathBW, '.tif', '_Edges_Canny.tif');
    imwrite(edgesCannyClean, edgeCannyPath);
    figure; imshow(edgesCannyClean);
    title('Canny Edge detection on binarised image');
    fprintf('Canny edge image exported: %s\n', edgeCannyPath);

    %% Extract coordinates of all edge points
    [edgeY, edgeX] = find(edgesCannyClean);

    %% Overlay traced outer edges and spray angle
    figure; imshow(BW_sauvola); hold on;
    distances_pix = round(distances_mm / pixelSize_mm);
    yDetect = yBottom + distances_pix;
    boundaryPointsLeft  = nan(length(yDetect), 2);
    boundaryPointsRight = nan(length(yDetect), 2);

    [rows, cols] = size(edgesCannyClean);
    for k = 1:length(yDetect)
        yline = yDetect(k) - yBottom + 1;
        plot([1, cols], [yline, yline], 'g--', 'LineWidth', 1.5);
        idx = find(edgeY == yline);
        if ~isempty(idx)
            boundaryPointsLeft(k,:)  = [min(edgeX(idx)), yline];
            boundaryPointsRight(k,:) = [max(edgeX(idx)), yline];
            plot(boundaryPointsLeft(k,1),  boundaryPointsLeft(k,2),  'bo', 'MarkerFaceColor','b', 'MarkerSize',6);
            plot(boundaryPointsRight(k,1), boundaryPointsRight(k,2), 'ro', 'MarkerFaceColor','r', 'MarkerSize',6);
        end
    end

    %% Compute spray angles
    validLeft  = find(~isnan(boundaryPointsLeft(:,1)));
    validRight = find(~isnan(boundaryPointsRight(:,1)));

    if length(validLeft) >= 2
        dxL = boundaryPointsLeft(validLeft(2),1) - boundaryPointsLeft(validLeft(1),1);
        dyL = boundaryPointsLeft(validLeft(2),2) - boundaryPointsLeft(validLeft(1),2);
        sprayAngleLeft = atan2d(abs(dxL), dyL);
        fprintf('Left edge angle: %.2f deg\n', sprayAngleLeft);
        plot(boundaryPointsLeft(validLeft,1), boundaryPointsLeft(validLeft,2), 'b-', 'LineWidth', 2);
    end

    if length(validRight) >= 2
        dxR = boundaryPointsRight(validRight(2),1) - boundaryPointsRight(validRight(1),1);
        dyR = boundaryPointsRight(validRight(2),2) - boundaryPointsRight(validRight(1),2);
        sprayAngleRight = atan2d(abs(dxR), dyR);
        fprintf('Right edge angle: %.2f deg\n', sprayAngleRight);
        plot(boundaryPointsRight(validRight,1), boundaryPointsRight(validRight,2), 'r-', 'LineWidth', 2);
    end

    if exist('sprayAngleLeft','var') && exist('sprayAngleRight','var')
        totalConeAngle = sprayAngleLeft + sprayAngleRight;
        fprintf('Total spray cone angle: %.2f deg\n', totalConeAngle);
    end

    title('Outer Spray Edges with Measurement Lines and Spray Angles');
    hold off;
end
%%
% 1. Raw background (cropped)
subplot(2,2,1);
imshow(bgMono(yBottom:end, :),[]);
title('Raw Background (Cropped)');

% 2. Averaged image with ROI rectangle (cropped)
subplot(2,2,2);
imshow(avgMonoCropped);
title('Averaged Image (Cropped)');
hold on;

% 3. Binary mask overlay (cropped)
subplot(2,2,3);
imshow(avgMonoCropped);
title('Binarised Mask (Cropped)');
hold on;

hMask = imshow(BW_sauvola);
set(hMask, 'AlphaData', 0.3);  % Adjust transparency

% 4. Spray edges (cropped)
subplot(2,2,4);
imshow(BW_sauvola, []); hold on;
visboundaries(edgesCannyClean, 'Color','r');
title('Spray with Canny Edge Overlay (Cropped)');





% Plot spray angle lines over averaged grayscale image
figure;
imshow(avgMonoCropped); hold on;

% Plot left boundary points
validLeft = find(~isnan(boundaryPointsLeft(:,1)));
if ~isempty(validLeft)
    plot(boundaryPointsLeft(validLeft,1), boundaryPointsLeft(validLeft,2), 'bo', ...
         'MarkerFaceColor','b','MarkerSize',6);
    plot(boundaryPointsLeft(validLeft,1), boundaryPointsLeft(validLeft,2), 'b-', ...
         'LineWidth',2);
end

% Plot right boundary points
validRight = find(~isnan(boundaryPointsRight(:,1)));
if ~isempty(validRight)
    plot(boundaryPointsRight(validRight,1), boundaryPointsRight(validRight,2), 'ro', ...
         'MarkerFaceColor','r','MarkerSize',6);
    plot(boundaryPointsRight(validRight,1), boundaryPointsRight(validRight,2), 'r-', ...
         'LineWidth',2);
end

% Title with total cone angle if both angles exist
if exist('sprayAngleLeft','var') && exist('sprayAngleRight','var')
    totalConeAngle = sprayAngleLeft + sprayAngleRight;
    title(sprintf('Spray Angle Overlay (Total Cone: %.2f deg)', totalConeAngle));
end

hold off;
