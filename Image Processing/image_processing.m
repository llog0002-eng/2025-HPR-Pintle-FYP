clear all; close all; clc;

%% --- User settings ---
processGif = true;
processSprayAngle = true;
allFrames = false;
monoChannel = 1;  % 1=Red, 2=Green, 3=Blue
folder = '\\ad.monash.edu\home\User090\ppin0001\Desktop\Pintle Data\19_9_25\';
filePath = append(folder,'1_1_3_1_60bar.mraw');
backgroundPath = '\\ad.monash.edu\home\User090\ppin0001\Desktop\Pintle Data\19_9_25\1_1_3_1_60bar.mraw';
outputPath16 = append(folder, 'avgMonoSubFrame.tif');
outputPathBW = append(folder,'avgMonoSubFrameBW.tif');
gifPath = append(folder,'bitshift.gif');

sigmaSmooth = 1;
backgroundFrame = 1;
maxFrames = 500;
bitshiftAmount = 2;
maxPixelValue = 4095;
pintleDiameter_mm = 25;
background_threshold = 0.009;

% Distances from the pintle tip for spray angle calculation
distances_mm = [0, 20];

rowMid = 50;
epsVal = 1e-6;
minPixelSize = 20;

% Camera info
frameWidth = 1024; frameHeight = 640; numChannels = 3; bitDepth = 12;

thresholdDarkFraction = 0.08;
startFrame = 1;
maxFrames = 500;

%% --- Load background ---
bgData = double(readmraw(backgroundPath, backgroundFrame));
bgMono = bgData(:,:,monoChannel);

startFrame = 1;             

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

    disp(['Preview GIF saved to: ', gifPath]);
end

% --- User-defined parameters ---
startFrame = 1;        % first frame to check
maxFrames = 400;         % last frame to check
frameJump = 1;           % jump for forward search (e.g., check every 5 frames)
numBgFrames = 200;        % number of frames to use for background variance
sigmaFactor = 2;         % threshold = mean + 3*std
maxPixelValue = 4095;    % 12-bit camera max value
monoChannel = 1;         % 1=Red, 2=Green, 3=Blue

exactFlowFrame = detectFlowFrames(filePath, monoChannel, startFrame, maxFrames, sigmaFactor, maxPixelValue);

disp(['Flow starts at frame: ', num2str(exactFlowFrame)]);


%% --- Display GIF (optional) ---
if processGif

    frame_start = exactFlowFrame-10;  % use detected flow start
    I = double(readmraw(filePath, frame_start));
    Ibright = uint16(I);
    Ibright = bitshift(Ibright, bitshiftAmount);
    Ibright(Ibright > maxPixelValue) = maxPixelValue;

    % Convert to grayscale for GIF
    Igray = uint8(255 * mat2gray(rgb2gray(Ibright)));

    h = imshow(Igray); % initial display
    imwrite(Igray, gifPath, 'gif', 'LoopCount', Inf, 'DelayTime', frameDelay); % first frame

    for i = frame_start+1 : frame_start + visualFrames - 1
        I = double(readmraw(filePath, i));

        % Bitshift all RGB channels
        Ibright = uint16(I);
        Ibright = bitshift(Ibright, bitshiftAmount);
        Ibright(Ibright > maxPixelValue) = maxPixelValue;

        % Convert to grayscale
        Igray = uint8(255 * mat2gray(rgb2gray(Ibright)));

        % Update display
        set(h, 'CData', Igray);
        drawnow;

        % Append frame to GIF
        imwrite(Igray, gifPath, 'gif', 'WriteMode', 'append', 'DelayTime', frameDelay);
    end

    disp(['GIF saved to: ', gifPath]);
end
%% --- Pintle detection (unchanged) ---
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
disp(['Detected pintle center at (X,Y) = (', num2str(centerX), ',', num2str(centerY), ')']);

colRange = startIdx:endIdx;
colData = any(bwPost(:,colRange),2);
yTop = find(colData,1,'first');
yBottom = find(colData,1,'last');
disp(['Pintle vertical extent from Y = ', num2str(yTop), ' to Y = ', num2str(yBottom)]);

pintleWidthPixels = endIdx - startIdx;
pixelSize_mm = pintleDiameter_mm / pintleWidthPixels;
disp(['Pixel size: ', num2str(pixelSize_mm), ' mm/pixel']);

%% --- Convert distances to pixels ---
distances_pix = round(distances_mm / pixelSize_mm);
yDetect = yBottom + distances_pix;

%% --- File info ---
fid = fopen(filePath,'r'); fseek(fid,0,'eof'); fileBytes = ftell(fid); fclose(fid);
bytesPerPixel = ceil(bitDepth/8);
bytesPerFrame = frameWidth * frameHeight * numChannels * bytesPerPixel;
numFrames = floor(fileBytes / bytesPerFrame);
numFrames = min(numFrames, maxFrames);
disp(['Total frames used: ', num2str(numFrames)]);

%% --- Frame averaging and background division ---
if processSprayAngle
    accumMono = zeros(frameHeight, frameWidth, 'double');

    % Accumulate over chunks
    for f = firstFlowFrame:chunkSize:numFrames
        endFrame = min(f + chunkSize - 1, numFrames);
        I_chunk = readmraw(filePath, [f endFrame]);

        % Background division for the mono channel
        monoChunk = double(I_chunk(:,:,monoChannel,:)) ./ (bgMono + epsVal);

        % Sum over 4th dimension (frames in chunk)
        accumMono = accumMono + sum(monoChunk, 4);
    end

    % Average over all frames
    avgMono = accumMono / numFrames;
%% --- Crop the averaged image like before ---
avgMonoCropped = avgMono(yBottom:end, :);
pct = 35; 
thresh = prctile(avgMonoCropped(:), pct);
BW_core = avgMonoCropped > thresh; 

window = [51 51];
k = 0.1;
BW_edges = sauvolaSingle(avgMonoCropped, window, k);
BW_edges = bwareaopen(BW_edges, minPixelSize);
BW_clean= BW_core & BW_edges;

figure;
imshow(BW_clean, []);
title('Combined threshold')


%% --- Save combined mask ---
imwrite(BW_clean, 'avgMonoHybridMask.jpg');
disp('Hybrid mask saved as avgMonoHybridMask.jpg');
    %% --- Edge detection ---
    edgesCanny = edge(BW_clean, 'Canny', [], 1);
    edgesCannyClean = bwareaopen(edgesCanny, 20);

    edgeCannyPath = strrep(outputPathBW, '.tif', '_Edges_Canny.tif');
    imwrite(edgesCannyClean, edgeCannyPath);
    disp(['Canny edge image exported: ', edgeCannyPath]);

    %% --- Spray angle calculation ---
    [edgeY, edgeX] = find(edgesCannyClean);
    boundaryPoints = zeros(2, 2);

    for kLine = 1:2
        yLine = yDetect(kLine);
        idx = find(edgeY == yLine);

        if isempty(idx)
            boundaryPoints(kLine, :) = [NaN NaN];
        else
            xIntersect = max(edgeX(idx));
            boundaryPoints(kLine, :) = [xIntersect, yLine];
        end
    end

    if all(~isnan(boundaryPoints(:)))
        vecLine = boundaryPoints(2, :) - boundaryPoints(1, :);
        sprayAngle = atan2d(vecLine(1), vecLine(2));
        disp(['Spray angle between ', num2str(distances_mm(1)), ' mm and ', ...
            num2str(distances_mm(2)), ' mm from pintle is ', num2str(sprayAngle), ' deg']);
    end
end


% Figure 1: Background, averaged, binary, and edges
figure('Name','Edge Detection Stages','NumberTitle','off');

% 1. Raw background
subplot(2,2,1);
imshow(mat2gray(bgMono));
title('Raw Background');

%% --- Subplot 2: Averaged image with ROI rectangle ---
subplot(2,2,2);
imshow(avgMono);  % show normalized averaged image
title('Averaged Image (Normalized by Background)');
hold on;

% Overlay ROI rectangle
yROI_start = yBottom + distances_pix(1);
yROI_end   = min(yBottom + distances_pix(2), frameHeight);
rectangle('Position', [1, yROI_start, frameWidth, yROI_end - yROI_start + 1], ...
    'EdgeColor', 'r', 'LineWidth', 2, 'LineStyle', '--');

%% --- Subplot 3: Binary mask overlay ---
subplot(2,2,3);
imshow(avgMono, []);  % show background-subtracted averaged image
title('Binarised mask')
hold on;

% Pad BW_clean back to full image size
BW_full = false(size(avgMono));
BW_full(yBottom:end, :) = BW_clean;  % place cropped mask at correct position

% Create AlphaData for transparency
alphaMask = zeros(size(avgMono));
alphaMask(yBottom:end, :) = 0.3;

% Overlay binary mask
hMask = imshow(BW_full);
set(hMask, 'AlphaData', alphaMask);

% 4. Spray edges
subplot(2,2,4);
imshow(BW_clean, []); hold on;
visboundaries(edgesCannyClean, 'Color','r');
title('Spray with Canny Edge Overlay');


%% --- Save images ---
if processSprayAngle
    % --- Save averaged image as JPG ---
    imwrite(uint8(255 * mat2gray(avgMono)), 'avgMono.jpg');
    disp('Averaged image exported as JPG: avgMono.jpg');

    imwrite(BW_clean, outputPathBW); % Fixed: previously used undefined avgMonoBW
    disp(['Binarized image exported: ', outputPathBW]);
end


function BW = sauvolaSingle(I, window, k)
% SAUVOLASINGLE - Single-stage Sauvola binarization
% I : grayscale image (normalized or raw)
% window : [M N] neighborhood (default [3 3])
% k : Sauvola sensitivity (default 0.34)

if nargin < 3
    k = 0.34;
end
if nargin < 2
    window = [3 3];
end

% Ensure input is 2D
if ~ismatrix(I)
    error('Input image must be a 2D array.');
end

I = double(I);

% Compute local mean and standard deviation using averagefilter
m  = averagefilter(I, window, 'replicate');
m2 = averagefilter(I.^2, window, 'replicate');
s  = sqrt(max(m2 - m.^2, 0));

% Sauvola threshold
R = max(s(:));
T = m .* (1 + k * (s./R - 1));

% Apply threshold
BW = I > T;
end

function image=averagefilter(image, varargin)
%AVERAGEFILTER 2-D mean filtering.
%   B = AVERAGEFILTER(A) performs mean filtering of two dimensional
%   matrix A with integral image method. Each output pixel contains
%   the mean value of the 3-by-3 neighborhood around the corresponding
%   pixel in the input image.
%
%   B = AVERAGEFILTER(A, [M N]) filters matrix A with M-by-N neighborhood.
%   M defines vertical window size and N defines horizontal window size.
%
%   B = AVERAGEFILTER(A, [M N], PADDING) filters matrix A with the
%   predefinned padding. By default the matrix is padded with zeros to
%   be compatible with IMFILTER. But then the borders may appear distorted.
%   To deal with border distortion the PADDING parameter can be either
%   set to a scalar or a string:
%       'circular'    Pads with circular repetition of elements.
%       'replicate'   Repeats border elements of matrix A.
%       'symmetric'   Pads array with mirror reflections of itself.
%
%   Comparison
%   ----------
%   There are different ways how to perform mean filtering in MATLAB.
%   An effective way for small neighborhoods is to use IMFILTER:
%
%       I = imread('eight.tif');
%       meanFilter = fspecial('average', [3 3]);
%       J = imfilter(I, meanFilter);
%       figure, imshow(I), figure, imshow(J)
%
%   However, IMFILTER slows down with the increasing size of the
%   neighborhood while AVERAGEFILTER processing time remains constant.
%   And once one of the neighborhood dimensions is over 21 pixels,
%   AVERAGEFILTER is faster. Anyway, both IMFILTER and AVERAGEFILTER give
%   the same results.
%
%   Remarks
%   -------
%   The output matrix type is the same as of the input matrix A.
%   If either dimesion of the neighborhood is even, the dimension is
%   rounded down to the closest odd value.
%
%   Example
%   -------
%       I = imread('eight.tif');
%       J = averagefilter(I, [3 3]);
%       figure, imshow(I), figure, imshow(J)
%
%   See also IMFILTER, FSPECIAL, PADARRAY.
%   Contributed by Jan Motl (jan@motl.us)
%   $Revision: 1.2 $  $Date: 2013/02/13 16:58:01 $
% Parameter checking.
numvarargs = length(varargin);
if numvarargs > 2
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 2 optional inputs');
end

optargs = {[3 3] 0};            % set defaults for optional inputs
optargs(1:numvarargs) = varargin;
[window, padding] = optargs{:}; % use memorable variable names
m = window(1);
n = window(2);
if ~mod(m,2) m = m-1; end       % check for even window sizes
if ~mod(n,2) n = n-1; end
if (ndims(image)~=2)            % check for color pictures
    display('The input image must be a two dimensional array.')
    display('Consider using rgb2gray or similar function.')
    return
end
% Initialization.
[rows columns] = size(image);   % size of the image
% Pad the image.
imageP  = padarray(image, [(m+1)/2 (n+1)/2], padding, 'pre');
imagePP = padarray(imageP, [(m-1)/2 (n-1)/2], padding, 'post');
% Always use double because uint8 would be too small.
imageD = double(imagePP);
% Matrix 't' is the sum of numbers on the left and above the current cell.
t = cumsum(cumsum(imageD),2);
% Calculate the mean values from the look up table 't'.
imageI = t(1+m:rows+m, 1+n:columns+n) + t(1:rows, 1:columns)...
    - t(1+m:rows+m, 1:columns) - t(1:rows, 1+n:columns+n);
% Now each pixel contains sum of the window. But we want the average value.
imageI = imageI/(m*n);
% Return matrix in the original type class.
image = cast(imageI, class(image));
end
