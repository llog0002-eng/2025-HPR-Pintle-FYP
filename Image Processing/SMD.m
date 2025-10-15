clear all; close all; clc;

%% --- User settings ---

% Set plots to light mode
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";

monoChannel = 1;  % 1=Red, 2=Green, 3=Blue

folder = 'A:\Uni\FYP\Droplet test data\';
datafolder = 'A:\Uni\FYP\Droplet test data\';
filePath = append(datafolder,'4_1_3_1_88bar.mraw');
pintlePath = append(datafolder,'4_1_3_1_88bar.mraw');
outfile = "SMDout.csv";

outputPath16 = append(folder, 'avgMonoSubFrame.tif');
outputPathBW = append(folder,'avgMonoSubFrameBW.tif');
gifPath = append(folder,'bitshift.gif');


frameDelay = 0.005; % delay between frames in seconds (for plotting)
visualFrames = 50; % number of frames to display in MATLAB plot (from start frame of flow)
chunkSize = 100; % allocation for efficiency 
percentileThresholdAng = 45; % threshold for detecting spray, spray angle
pixelSize_mm = 0.0996;

sigmaSmooth = 1;
startFrame = 4;
sigmaFactor = 3;
backgroundFrame = 1; % frame to detect background
frame_start = 200; % frame to detect initial flow through pintle
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
frameNoSMDstart = 650; % Frame number to start analysing droplets
frameNoSMDstep = 10; % Frames to step between droplet analysis
frameNoSMDend = 750; % Frame number to stop analysing droplets
maxSMD = 2; % Max SMD to plot for histogram, mm
percentileThresholdSMD = 20; % threshold for detecting spray, SMD
connectivity = 4;

% Cropping inputs
% Crop to region of interest (ROI)
cropHeight = 400; 
cropWidth  = 400;
yCenter = 400;
xCenter = 500;

% Camera settings (only used for determine info in command window not used
% in script)

frameWidth = 1024; frameHeight = 640; numChannels = 3; bitDepth = 12;

%% Pintle Post Detection

%exactFlowFrame = detectFlowFrames(filePath, monoChannel, startFrame, maxFrames, sigmaFactor, maxPixelValue);

% --- Read background frame ---
bgData = double(readmraw(pintlePath, 1));
pintleMono = bgData(:,:,monoChannel);

bwPost = pintleMono < background_threshold*max(pintleMono(:)); % thresholds maximum values in background image
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
%pixelSize_mm = pintleDiameter_mm / pintleWidthPixels;
disp(['Pixel size: ', num2str(pixelSize_mm), ' mm/pixel']);

%% File info

fid = fopen(filePath,'r'); fseek(fid,0,'eof'); fileBytes = ftell(fid); fclose(fid);
bytesPerPixel = ceil(bitDepth/8);
bytesPerFrame = frameWidth * frameHeight * numChannels * bytesPerPixel;
numFrames = floor(fileBytes / bytesPerFrame);
numFrames = min(numFrames, maxFrames);
disp(['Total frames used: ', num2str(numFrames)]);

%% Thresholding

% Preallocate results:
frames = zeros([cropWidth,cropHeight]);
I_enhanceds = zeros([cropWidth,cropHeight]);
BWs = zeros([cropWidth,cropHeight]);

n = (frameNoSMDend - frameNoSMDstart)/frameNoSMDstep + 1;
dropletDiameters = NaN([1000,n]); % Preallocate NaNs to not use preallocation values
dropletAreas = NaN([1000,n]); % Preallocate NaNs to not use preallocation values

% For each frame being analysed:
for i=frameNoSMDstart:frameNoSMDstep:frameNoSMDend
    j = (i-frameNoSMDstart)/frameNoSMDstep+1; % Frame index

    % Read frames
    bg = readmraw(filePath, backgroundFrame);
    bgMono = double(bg(:,:,monoChannel));
    bgMonoNorm = mat2gray(bgMono);
    bg = bgMonoNorm;
    
    frame = readmraw(filePath, i);
    frameMono = double(frame(:,:,monoChannel));
    frameMonoNorm = mat2gray(frameMono);
    
    % Cropping
    y1 = max(1, round(yCenter - cropHeight/2));
    y2 = min(size(frameMonoNorm,1), round(yCenter + cropHeight/2));
    x1 = max(1, round(xCenter - cropWidth/2));
    x2 = min(size(frameMonoNorm,2), round(xCenter + cropWidth/2));
    
    frameMonoNorm = frameMonoNorm(y1:y2, x1:x2);
    bg = bg(y1:y2,x1:x2);
    
    % Subtract background
    frameMonoNorm = frameMonoNorm ./ bg;
    
    %% Matlab inbuilt processing
    frame = frameMonoNorm;
    I = frame;                % Convert to grayscale
    I = im2double(I);                  % Convert to double
    
    I_gray = im2double(im2gray(frame));        % convert
    % Smooth and enhance
    I_gray = imguidedfilter(I_gray, 'DegreeOfSmoothing', 0.005); 
    I_gray = adapthisteq(I_gray);  % optional CLAHE
    
    % Compute enhancement without scaling
    I_enhanced = 1 - I_gray; % Droplets are white
    I_enhanced(I_gray < 0) = 0; % clip negatives
    
    % Clip extreme values to avoid saturation
    I_enhanced(I_enhanced > 1) = 1;
    
    % Adaptive threshold
    T = adaptthresh(I_enhanced, 0.4, 'ForegroundPolarity','bright', 'NeighborhoodSize', 11);
    BW = imbinarize(I_enhanced, T);
    
    % Clean mask
    BW = bwareaopen(BW, 2, connectivity);  % remove small noise
    BW = imfill(BW, 'holes');
    
    % Label droplets
    CC = bwconncomp(BW,connectivity);
    
    % Measure droplet properties
    stats = regionprops(CC, 'Area', 'Centroid', 'BoundingBox');
    
    % Compute droplet diameters
    dropletArea = [stats.Area];
    dropletDiameters(1:length(dropletArea),j) = 2*sqrt(dropletArea/pi);  % Diameter in pixels
    dropletAreas(1:length(dropletArea),j) = dropletArea; % Area in pixels
    
    fprintf('Detected %d droplets.\n', length(stats));
    %fprintf('Estimated Mean Droplet Diameter = %.2f pixels.\n', mean(dropletDiameters(1:length(dropletAreas),j)));
end

D20s = pixelSize_mm*sqrt(4*dropletAreas(:)/pi); % surface diameters, D20
D20 = mean(D20s,"omitnan");

%% Plotting

figure('Name','Droplet Analysis','NumberTitle','off');
t = tiledlayout(2,2,'Padding','compact','TileSpacing','compact'); % 4 images, compact spacing

% Original (cropped) frame
nexttile;
imshow(frame); 
title('Original Cropped Frame');

% Enhanced Image
nexttile;
imshow(I_enhanced); 
title('Enhanced Image');

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

% % Area in pixels
% figure('Name','Droplet Size Distribution','NumberTitle','off');
% 
% histogram(dropletAreas(:), "BinMethod", "integers", "BinLimits", [0, 50]); % histogram with pixel-sized bins limited to 50
% xlabel('Droplet Area (pixels)');
% ylabel('Number of Droplets');
% title('Droplet Size Distribution');
% grid on;

% Diameter in mm
figure('Name','Droplet Size Distribution','NumberTitle','off');

edges = linspace(0,2,25);
histogram(D20s, edges); % histogram with pixel-sized bins limited to 50
xlabel('Droplet Surface Diameter (mm)');
ylabel('Number of Droplets');
title('Droplet Size Distribution');
grid on;

% Just droplet centroids on original image

figure;
imshow(frame); 
hold on;
for k = 1:length(stats)
    plot(stats(k).Centroid(1), stats(k).Centroid(2), 'r*', 'MarkerSize', 5);
end
hold off;
title('Droplet Centroids on Original Image');

%% Text output

fprintf("D20: %.2f mm\n",D20)

%% File Output

hpt = 1.02; % Height of pintle, mm, to write to first column
writematrix([hpt D20s'],outfile,'Delimiter',',','WriteMode','append');
