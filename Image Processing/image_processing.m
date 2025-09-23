clear all; close all; clc;

%% --- User settings ---
% Toggle features
processGif = false;
processSprayAngle = false;
processSMD = true;

monoChannel = 1;  % 1=Red, 2=Green, 3=Blue
folder = 'A:\OneDrive - Monash University\Uni\HPR\FYP\Testing data\Sample\';
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
pintley = 200; % For cropping out pintle

% SMD
frameNoSMD = 400; % Frame number to analyse SMD
maxSMD = 2; % Max SMD to plot for histogram, mm
percentileThresholdSMD = 20; % threshold for detecting spray, SMD

% Camera settings (only used for determine info in command window not used
% in script)

frameWidth = 1024; frameHeight = 640; numChannels = 3; bitDepth = 12;

%% --- Display bitshifted raw frames and save as GIF ---

if processGif == true
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


%% Pintle detection on background 
% --- Read background frame ---
bgData = double(readmraw(filePath, backgroundFrame));
bgMono = bgData(:,:,monoChannel);
figure('Name','Raw Background Image','NumberTitle','off');
imshow(mat2gray(bgMono));
title('Raw Background Image');

bwPost = bgMono < background_threshold*max(bgMono(:)); % thresholds maximum values in background image
bwPost = imfill(bwPost,'holes'); 
bwPost = bwareaopen(bwPost,50); 

rowData = bwPost(rowMid,:);
figure('Name','Row Data Visualization','NumberTitle','off');
plot(rowData,'LineWidth',2); ylim([-0.1 1.1]);
xlabel('X [pixels]'); ylabel('Pixel Value (0=black, 1=white)');
title(['Row Data at Y = ', num2str(rowMid)]);

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

figure('Name','Thresholded Background (Pintle Detection)','NumberTitle','off');
imshow(bwPost); hold on;
rectangle('Position', pintleBox, 'EdgeColor', 'c', 'LineWidth', 2);
yVals = 1:size(bgMono,1);
plot(centerX*ones(size(yVals)), yVals, 'b--', 'LineWidth', 2);
text(centerX+5, centerY, 'Pintle Center', 'Color', 'y', 'FontSize', 12);
title('Thresholded Background with Pintle Center');
axis on; axis image; xlabel('X [pixels]'); ylabel('Y [pixels]'); grid on;

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

%% Frame averaging and division by background (one channel)'
if processSprayAngle == true
    
    accumMono = zeros(frameHeight, frameWidth, 'double');
    
    for startFrame = frame_start:chunkSize:numFrames
        endFrame = min(startFrame+chunkSize-1, numFrames);
        I_chunk = readmraw(filePath, [startFrame endFrame]);
        
        % Divide by background (reads one channel)
        monoChunk = double(I_chunk(:,:,monoChannel,:)) ./ (bgMono + epsVal);
      
        accumMono = accumMono + sum(monoChunk, 4);
    end
    
    avgMono = accumMono / numFrames; % average frame
    avgMonoNorm = mat2gray(avgMono); 
    
    %% Thresholding using percentile-based method
    % Compute threshold based on percentile of pixel intensities
    threshVal = prctile(avgMonoNorm(:), percentileThresholdAng);
    
    % Binarize image
    avgMonoBW = avgMonoNorm < threshVal;
    
    % Morphological cleaning
    avgMonoBW = bwareaopen(avgMonoBW, 50);         
    avgMonoBW = imclose(avgMonoBW, strel('disk',3)); 
    
    %% Edge detection
    sigma = 1;        
    threshold = [];   
    
    edgesCanny = edge(avgMonoBW, 'Canny', threshold, sigma);
    edgesCannyClean = bwareaopen(edgesCanny, 20);  
    
    edgeCannyPath = strrep(outputPathBW, '.tif', '_Edges_Canny.tif');
    imwrite(edgesCannyClean, edgeCannyPath);
    disp(['Canny edge image exported: ', edgeCannyPath]);
    
    %% Spray angle
    pintleCenter = [centerX, centerY];
    [edgeY, edgeX] = find(edgesCannyClean);
    boundaryPoints = zeros(2,2);
    
    for k = 1:2
        yLine = yDetect(k);
        idx = find(edgeY == yLine);
        if isempty(idx)
            disp(['No edge pixels found at y = ' num2str(yLine)]);
            boundaryPoints(k,:) = [NaN NaN];
        else
            xIntersect = max(edgeX(idx)); % right side
            boundaryPoints(k,:) = [xIntersect, yLine];
        end
    end
    
    if all(~isnan(boundaryPoints(:)))
        vecLine = boundaryPoints(2,:) - boundaryPoints(1,:);
        sprayAngle = atan2d(vecLine(1), vecLine(2));
        disp(['Spray angle between y = ' num2str(yVals(1)) ...
              ' and y = ' num2str(yVals(2)) ' is ' ...
              num2str(sprayAngle) ' deg']);
    
        % --- Visualization ---
        figure('Name','Spray Line from Two Intersections','NumberTitle','off');
        imshow(avgMonoBW); hold on;
        yline(centerY,'b','DisplayName','Line used to detect pintle')
        plot(boundaryPoints(:,1), boundaryPoints(:,2), 'rx', 'MarkerSize', 10, 'LineWidth', 2,'DisplayName','Intersection Points');
        line(boundaryPoints(:,1), boundaryPoints(:,2), 'Color','r','LineWidth',2,'DisplayName','Line between points');
        line([pintleCenter(1), pintleCenter(1)], ...
             [pintleCenter(2), size(avgMono,1)], 'Color','m','LineStyle','--','LineWidth',1.5,'DisplayName','Pintle Centerline');
       yline(yDetect,'g-','DisplayName','Pintle Angle Detection Points')
        title(sprintf('Spray Angle = %.2f°', sprayAngle));
        legend()
    end
end

%% SMD

if processSMD == true
    
    frame = readmraw(filePath, frameNoSMD);
    frameMono = double(frame(:,:,monoChannel));
    frameMonoNorm = mat2gray(frameMono);

    % Subtract background
    frameMonoNorm = frameMonoNorm ./ mat2gray(bgMono);

    figure();
    imshow(frameMonoNorm);

    % Crop out pintle
    frameMonoNorm = frameMonoNorm(pintley:end,:);

    % Binarise image
    % Compute threshold based on percentile of pixel intensities
    threshVal = prctile(frameMonoNorm(:), percentileThresholdSMD);
    frameBW = frameMonoNorm < threshVal;

    figure;
    imshow(frameBW);

    % Clustering
    % Initialize label matrix
    [labelMatrix, numLabels] = bwlabel(frameBW, 8); % 8-connectivity
    % Display the labeled image with a custom colormap (zero as black)
    figure;
    imshow(labelMatrix, []);
    title(['Labeled Image with ', num2str(numLabels), ' Regions']);
    colormap([0 0 0; hsv(numLabels)]); % Apply custom colormap with zero as black
    colorbar; % Add colorbar for reference

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

    figure;
    imshow(clusterImage, []);
    colormap([0 0 0; parula(numLabels)]);
    SMDcolorbar=colorbar; % Add colorbar for reference
    SMDcolorbar.Label.String = "SMD, mm";
    
    figure;
    binsSMD = maxSMD/(4/pi*pixelSize_mm);
    SMDedges = linspace(0,maxSMD,binsSMD);
    histogram(cellEffDia,"BinEdges",SMDedges)
    title("Particle Diameter Distribution")
    xlabel("SMD, mm")
    ylabel("Occurrences")
    
    % Overlay cluster image on the background-subtracted frame
    figure;
    imshow(frameMonoNorm, []); hold on;
    % Create a binary mask for the cluster image
    clusterMask = clusterImage > 0; 
    % Overlay red outlines for the clusters
    visboundaries(clusterMask, 'Color', 'r', 'LineWidth', 1);
    title('Cluster Image Overlay on Background-Subtracted Frame');
end

%% Plotting
if processSprayAngle == true
    figure('Name','Edge Detection Stages','NumberTitle','off');
    
    subplot(2,2,1);
    imshow(mat2gray(bgMono)); 
    title('Raw Background');
    
    subplot(2,2,2)
    imshow(avgMono)
    title('Averaged image (normalised by background)')
    subplot(2,2,3);
    imshow(avgMonoBW, []);  
    title('Averaged image (normalised by background) overlayed with threshold image'); 
    hold on;
    
    % Overlay binary mask in red with transparency
    h = imshow(avgMono);
    set(h, 'AlphaData', 0.7);  % 0 = transparent, 1 = opaque
    colormap(gca, 'gray');
    
    
    subplot(2,2,4);
    imshow(avgMono, []); hold on;  
    visboundaries(edgesCannyClean, 'Color','r'); 
    title('Spray with Canny Edge Overlay');
    
    figure('Name','Averaged background-subtracted image','NumberTitle','off');
    imshow(avgMonoBW, []); colormap('gray'); colorbar;
    title('Averaged background-subtracted (normalized 0–1)');
end

%% Saving images

if processSprayAngle == true
    imwrite(uint16(avgMonoNorm * maxPixelValue), outputPath16, 'tif');
    disp(['16-bit averaged frame exported: ', outputPath16]);
    
    imwrite(avgMonoBW, outputPathBW);
    disp(['Binarized image exported: ', outputPathBW]);
end
