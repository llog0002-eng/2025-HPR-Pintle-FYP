clear all; close all; clc;

%% --- User settings ---

monoChannel = 1;  % 1=Red, 2=Green, 3=Blue
folder = 'C:\Users\Luke''s Laptop\OneDrive - Monash University\Uni\HPR\FYP\Testing data\Sample\';
filePath = append(folder,'1_1_4_4_80bar.mraw');
backgroundPath = append(folder,'background.mraw');
outputPath16 = append(folder, 'avgMonoSubFrame.tif');
outputPathBW = append(folder,'avgMonoSubFrameBW.tif');
gifPath = append(folder,'bitshift.gif');
frameDelay = 0.005; % delay between frames in seconds (for plotting)
visualFrames = 50; % number of frames to display in MATLAB plot (from start frame of flow)
chunkSize = 100; % allocation for efficiency 
percentileThresholdAng = 45; % threshold for detecting spray, spray angle

sigmaSmooth = 1;
backgroundFrame = 1; % frame to detect background
frame_start = 100; % frame to detect initial flow through pintle
maxFrames = 200; % maximum number of frames to sample (from initial flow)
bitshiftAmount = 2; % used for visualisation of the pintle in GIF
maxPixelValue = 4095; % 12-bit max
pintleDiameter_mm = 25;  
background_threshold = 0.009; % threshold for detecting pintle (very low given pintle is dark)
yDetect = [150, 200];   % pintle spray angle measured two lines at these coordinates y = 0 is the top of the frame 
rowMid = 50; % y-coordinate for detecting pintle center
epsVal = 1e-6; % Background division epsilon value
pintley = 450; % For cropping out pintle

% SMD
frameNoSMD = 400; % Frame number to analyse SMD
maxSMD = 2; % Max SMD to plot for histogram, mm
percentileThresholdSMD = 20; % threshold for detecting spray, SMD
connectivity = 4;

% Cropping inputs
% Crop to region of interest (ROI)
cropHeight = 100; 
cropWidth  = 100;

yCenter = 500; % adjust as needed (towards bottom of 191 px image)
xCenter = 512; % roughly center of 1024 px wide image

% Camera settings (only used for determine info in command window not used
% in script)

frameWidth = 1024; frameHeight = 640; numChannels = 3; bitDepth = 12;

%% Background

% --- Read background frame ---
bgData = double(readmraw(filePath, backgroundFrame));
bgMono = bgData(:,:,monoChannel);

bwPost = bgMono < background_threshold*max(bgMono(:)); % thresholds maximum values in background image
bwPost = imfill(bwPost,'holes'); 
bwPost = bwareaopen(bwPost,50); 

rowData = bwPost(rowMid,:);

transitions = diff(rowData);
startIdx = find(transitions == 1, 1, 'first');
endIdx   = find(transitions == -1, 1, 'last');

if ~isempty(startIdx) && ~isempty(endIdx)
    centerX = round((endIdx - startIdx)/2 + startIdx);
    centerY = rowMid;
    pintleBox = [startIdx, rowMid-10, endIdx-startIdx, 20];
else
    centerX = size(bwPost,2)/2;
    centerY = rowMid;
    pintleBox = [centerX-1, centerY-1, 2, 2];
end
disp(['Detected pintle center at (X, Y) = (', num2str(centerX), ', ', num2str(centerY), ')']);



pintleWidthPixels = endIdx - startIdx;  % width of pintle in pixels
pixelSize_mm = pintleDiameter_mm / pintleWidthPixels;
disp(['Pixel size: ', num2str(pixelSize_mm), ' mm/pixel']);

%% File info

fid = fopen(filePath,'r'); fseek(fid,0,'eof'); fileBytes = ftell(fid); fclose(fid);
bytesPerPixel = ceil(bitDepth/8);
bytesPerFrame = frameWidth * frameHeight * numChannels * bytesPerPixel;
numFrames = floor(fileBytes / bytesPerFrame);
numFrames = min(numFrames, maxFrames);
disp(['Total frames used: ', num2str(numFrames)]);

%% Thresholding

frame = readmraw(filePath, frameNoSMD);
frameMono = double(frame(:,:,monoChannel));
frameMonoNorm = mat2gray(frameMono);

% Subtract background
frameMonoNorm = frameMonoNorm ./ mat2gray(bgMono);

% Cropping
y1 = max(1, round(yCenter - cropHeight/2));
y2 = min(size(frameMonoNorm,1), round(yCenter + cropHeight/2));
x1 = max(1, round(xCenter - cropWidth/2));
x2 = min(size(frameMonoNorm,2), round(xCenter + cropWidth/2));

frameMonoNorm = frameMonoNorm(y1:y2, x1:x2);
bg = bgMono(y1:y2,x1:x2);

% Binarise image
% Compute threshold based on percentile of pixel intensities
threshVal = prctile(frameMonoNorm(:), percentileThresholdSMD);
frameBW = frameMonoNorm < threshVal;

% Clustering
% Initialize label matrix
[labelMatrix, numLabels] = bwlabel(frameBW, connectivity); % 4-connectivity

% Calculate SMD
cellCounts = zeros(numLabels, 1); % Initialize array to hold cell counts
cellEffDia = zeros(numLabels, 1); % Initialize array to hold effective cell diameters
for i = 1:numLabels
    cellCounts(i) = sum(labelMatrix(:) == i); % Count number of pixels in each cluster
    cellEffDia(i) = sqrt(4*cellCounts(i)/pi)*pixelSize_mm; % Convert to effective cell dia (eq to circle), in mm
end

% Create an image of the clusters colored by effective cell diameter
clusterImage = zeros(size(frameBW, 1), size(frameBW, 2)); % Initialize image
for i = 1:numLabels
    if cellEffDia(i) < maxSMD
        clusterImage(labelMatrix==i) = cellEffDia(i);
    else
        clusterImage(labelMatrix==i) = 0;
    end
end

%% Plotting

% % Display background image
% figure('Name','Raw Background Image','NumberTitle','off');
% imshow(mat2gray(bgMono));
% title('Raw Background Image');

% % Display row data
% figure('Name','Row Data Visualization','NumberTitle','off');
% plot(rowData,'LineWidth',2); ylim([-0.1 1.1]);
% xlabel('X [pixels]'); ylabel('Pixel Value (0=black, 1=white)');
% title(['Row Data at Y = ', num2str(rowMid)]);

% % Display pintle detection:
% figure('Name','Thresholded Background (Pintle Detection)','NumberTitle','off');
% imshow(bwPost); hold on;
% rectangle('Position', pintleBox, 'EdgeColor', 'c', 'LineWidth', 2);
% yVals = 1:size(bgMono,1);
% plot(centerX*ones(size(yVals)), yVals, 'b--', 'LineWidth', 2);
% text(centerX+5, centerY, 'Pintle Center', 'Color', 'y', 'FontSize', 12);
% title('Thresholded Background with Pintle Center');
% axis on; axis image; xlabel('X [pixels]'); ylabel('Y [pixels]'); grid on;

% % Display background subtracted frame:
% figure();
% imshow(frameMonoNorm);

% % Display binarised frame:
% figure;
% imshow(frameBW);

% % Display the labeled image with a custom colormap (zero as black)
% figure;
% imshow(labelMatrix, []);
% title(['Labeled Image with ', num2str(numLabels), ' Regions']);
% colormap([0 0 0; hsv(numLabels)]); % Apply custom colormap with zero as black
% colorbar; % Add colorbar for reference

% figure;
% imshow(clusterImage, []);
% colormap([0 0 0; parula(numLabels)]);
% SMDcolorbar=colorbar; % Add colorbar for reference
% SMDcolorbar.Label.String = "SMD, mm";

% figure;
% binsSMD = maxSMD/(4/pi*pixelSize_mm);
% SMDedges = linspace(0,maxSMD,binsSMD);
% histogram(cellEffDia,"BinEdges",SMDedges)
% title("Particle Diameter Distribution")
% xlabel("SMD, mm")
% ylabel("Occurrences")

% % Overlay cluster image on the background-subtracted frame
% figure;
% imshow(frameMonoNorm, []); hold on;
% % Create a binary mask for the cluster image
% clusterMask = clusterImage > 0; 
% % Overlay red outlines for the clusters
% visboundaries(clusterMask, 'Color', 'r', 'LineWidth', 1);
% title('Cluster Image Overlay on Background-Subtracted Frame');

%% Computer Vision
% Read an image

%% Load image
frame = frameMonoNorm;
I = frame;                % Convert to grayscale
I = im2double(I);                  % Convert to double

I_gray = im2double(im2gray(frame));        % convert
% Smooth and enhance
I_gray = imguidedfilter(I_gray, 'DegreeOfSmoothing', 0.01); 
I_gray = adapthisteq(I_gray);  % optional CLAHE

% Compute enhancement without scaling
I_enhanced = bg - I_gray;       % droplets are positive
I_enhanced(I_enhanced < 0) = 0; % clip negatives
I_enhanced = I_enhanced * 5;    % boost contrast without scaling to 1

% Clip extreme values to avoid saturation
I_enhanced(I_enhanced > 1) = 1;

% Adaptive threshold
T = adaptthresh(I_enhanced, 0.5, 'ForegroundPolarity','bright', 'NeighborhoodSize', 35);
BW = imbinarize(I_enhanced, T);

% Clean mask
BW = bwareaopen(BW, 2);  % remove small noise
BW = imfill(BW, 'holes');

% Step 2: Edge detection using Canny
%BW = edge(I_smooth, 'Canny', [0.05 0.2]);

% Step 3: Morphological cleanup
BW_clean = bwareaopen(BW, 2);   % Remove small objects
BW_clean = imfill(BW_clean, 'holes');

% Step 4: Label droplets
CC = bwconncomp(BW_clean,connectivity);

% Step 5: Measure droplet properties
stats = regionprops(CC, 'Area', 'Centroid', 'BoundingBox');

% Step 7: Compute droplet diameters
dropletAreas = [stats.Area];
dropletDiameters = 2*sqrt(dropletAreas/pi);  % diameter in pixels

% Step 8: Compute Sauter Mean Diameter (D32)
SMDs = sum(dropletDiameters.^3) / sum(dropletDiameters.^2);

fprintf('Detected %d droplets.\n', length(stats));
fprintf('Mean droplet area = %.2f mm2.\n', mean(dropletAreas)*pixelSize_mm^2);
fprintf('Estimated Sauter Mean Diameter (SMD) = %.2f mm.\n', SMDs*pixelSize_mm);

%% Plotting

figure('Name','Droplet Analysis','NumberTitle','off');
t = tiledlayout(3,2,'Padding','compact','TileSpacing','compact'); % 4 images, compact spacing

% Original (cropped) frame
nexttile;
imshow(frame); 
title('Original Cropped Frame');

% Enhanced Image
nexttile;
imshow(I_enhanced); 
title('Enhanced Image');

% Threshold
nexttile;
imshow(T);
title("Threshold")

% Binary Mask
nexttile;
imshow(BW);
title('Binary Mask');

% Original Frame with Centroids
nexttile;
imshow(frame); 
hold on;
for k = 1:length(stats)
    plot(stats(k).Centroid(1), stats(k).Centroid(2), 'r*', 'MarkerSize', 5);
end
hold off;
title('Droplet Centroids on Original Image');

sgtitle('Droplet Analysis'); % shared title

% Histogram

figure('Name','Droplet Size Distribution','NumberTitle','off');

% Define bin edges based on pixel size
maxDiameter_mm = max(dropletDiameters)*pixelSize_mm;   % maximum droplet diameter in mm
binEdges = 0:pixelSize_mm:(maxDiameter_mm + pixelSize_mm); % bins of size = 1 pixel

histogram(dropletDiameters*pixelSize_mm, binEdges); % histogram with pixel-sized bins
xlabel('Droplet Diameter (mm)');
ylabel('Number of Droplets');
title('Droplet Size Distribution');
grid on;