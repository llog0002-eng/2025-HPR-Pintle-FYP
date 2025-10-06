clear all; close all; clc;

%% User settings
allFrames = true;
monoChannel = 2;  % 1=Red, 2=Green, 3=Blue
folder = '\\ad.monash.edu\home\User090\ppin0001\Desktop\Pintle Data\19_9_25\';
backgroundPath = fullfile(folder,'1_1_3_1_60bar.mraw'); % static background

% Specify output directory
outputDir = '\\ad.monash.edu\home\User090\ppin0001\Desktop\BatchDataTest';
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

% Excel output path
outputExcel = fullfile(outputDir,'SprayAnalysisResults.xlsx');
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

distances_mm = [5, 15]; % spray angle distances
sigmaFactor = 0.008;
windowSize = 101;
kValue = 0.3;
rowMid = 50;
epsVal = 1e-6;
minPixelSize = 60;
maxFiles = 3; % maximum number of .mraw files to process
frameWidth = 1024; frameHeight = 640; numChannels = 3; bitDepth = 12;

%% Load background
bgData = double(readmraw(backgroundPath, backgroundFrame));
bgMono = bgData(:,:,monoChannel);

%% Find all .mraw files in folder
mrawFiles = dir(fullfile(folder,'*.mraw'));
allResults = {}; % to store results

numFiles = min(length(mrawFiles), maxFiles);

for fIdx = 1:numFiles
    fileName = mrawFiles(fIdx).name;
    filePath = fullfile(folder,fileName);
    fprintf('Processing %s ...\n', fileName);

    baseName = erase(fileName,'.mraw'); % define before use
    %% --- Check if file can be read ---
    try
        % Attempt to read first frame to confirm header exists
        testFrame = readmraw(filePath, 1); %#ok<NASGU>
    catch ME
        if contains(ME.message, 'Could not locate .CIH or .CIHX header')
            warning('Skipping file %s: no .CIH/.CIHX header found.', fileName);
            continue; % skip to next file
        else
            rethrow(ME); % any other error, stop
        end
    end

    %% File info
    fid = fopen(filePath,'r'); fseek(fid,0,'eof'); fileBytes = ftell(fid); fclose(fid);
    bytesPerPixel = ceil(bitDepth/8);
    bytesPerFrame = frameWidth*frameHeight*numChannels*bytesPerPixel;
    numFrames = floor(fileBytes/bytesPerFrame);
    numFrames = min(numFrames, maxFrames);

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
    else
        centerX = size(bwPost,2)/2;
    end
    colRange = startIdx:endIdx;
    colData = any(bwPost(:,colRange),2);
    yBottom = find(colData,1,'last');

    pintleWidthPixels = endIdx - startIdx;
    pixelSize_mm = pintleDiameter_mm / pintleWidthPixels;

    %% Convert distances to pixels
    distances_pix = round(distances_mm / pixelSize_mm);
    yDetect = yBottom + distances_pix;

    %% --- Detect flow frames ---
    [exactFlowFrame, frameList, changeList] = detectFlowFrames(filePath, monoChannel, 1, maxFrames, sigmaFactor, maxPixelValue);




    %% --- Create GIF if requested ---
    gifPath = fullfile(outputDir, [baseName '_AllFrames.gif']);
    if allFrames
        startFrame = exactFlowFrame; % use detected flow frame
        I = double(readmraw(filePath, startFrame));
        Ibright = uint16(I);
        Ibright = bitshift(Ibright, bitshiftAmount);
        Ibright(Ibright > maxPixelValue) = maxPixelValue;
        Igray = uint8(255 * mat2gray(rgb2gray(Ibright)));
        imwrite(Igray, gifPath, 'gif', 'LoopCount', Inf, 'DelayTime', frameDelay);

        for i = startFrame+1 : min(startFrame + visualFrames - 1, numFrames)
            I = double(readmraw(filePath, i));
            Ibright = uint16(I);
            Ibright = bitshift(Ibright, bitshiftAmount);
            Ibright(Ibright > maxPixelValue) = maxPixelValue;
            Igray = uint8(255 * mat2gray(rgb2gray(Ibright)));
            imwrite(Igray, gifPath, 'gif', 'WriteMode', 'append', 'DelayTime', frameDelay);
        end
        fprintf('GIF saved to: %s\n', gifPath);
    end

    %% --- Frame averaging and background division ---
    accumMono = zeros(frameHeight, frameWidth, 'double');
    for cf = exactFlowFrame:chunkSize:numFrames
        endFrame = min(cf+chunkSize-1, numFrames);
        I_chunk = readmraw(filePath, cf:endFrame);
        monoChunk = double(I_chunk(:,:,monoChannel,:)) ./ bgMono;
        accumMono = accumMono + sum(monoChunk,4);
    end
    avgMono = accumMono / numFrames; % keep numFrames as divisor
    avgMonoCropped = imadjust(avgMono(yBottom:end,:));

    %% --- Sauvola mask ---
    BW_sauvola = sauvolaSingle(avgMonoCropped,[windowSize windowSize], kValue);
    BW_inv = ~BW_sauvola;
    BW_sauvola_inv = bwareaopen(BW_inv,5,4);
    BW_sauvola = ~BW_sauvola_inv;

    %% --- Edge detection ---
    edgesCanny = edge(BW_sauvola,'Canny',[],3);
    edgesCannyClean = bwareaopen(edgesCanny, minPixelSize);
    edgesCannyClean = imclose(edgesCannyClean, strel('line',5,0));

    %% --- Overlay traced edges ---
    [edgeY, edgeX] = find(edgesCannyClean);
    boundaryPointsLeft  = nan(length(yDetect),2);
    boundaryPointsRight = nan(length(yDetect),2);
    for k = 1:length(yDetect)
        yline = yDetect(k)-yBottom+1;
        idx = find(edgeY == yline);
        if ~isempty(idx)
            boundaryPointsLeft(k,:) = [min(edgeX(idx)), yline];
            boundaryPointsRight(k,:) = [max(edgeX(idx)), yline];
        end
    end

    validLeft = find(~isnan(boundaryPointsLeft(:,1)));
    validRight = find(~isnan(boundaryPointsRight(:,1)));

    if length(validLeft)>=2
        dxL = boundaryPointsLeft(validLeft(2),1)-boundaryPointsLeft(validLeft(1),1);
        dyL = boundaryPointsLeft(validLeft(2),2)-boundaryPointsLeft(validLeft(1),2);
        sprayAngleLeft = atan2d(abs(dxL), dyL);
    else
        sprayAngleLeft = NaN;
    end
    if length(validRight)>=2
        dxR = boundaryPointsRight(validRight(2),1)-boundaryPointsRight(validRight(1),1);
        dyR = boundaryPointsRight(validRight(2),2)-boundaryPointsRight(validRight(1),2);
        sprayAngleRight = atan2d(abs(dxR), dyR);
    else
        sprayAngleRight = NaN;
    end
    totalConeAngle = sprayAngleLeft + sprayAngleRight;

    %% --- Save images ---
    binarisationPath = fullfile(outputDir,[baseName '_SauvolaBinarisation.png']);
    edgePath = fullfile(outputDir,[baseName '_CannyEdges.png']);
    avgPath = fullfile(outputDir,[baseName '_AveragedImage.png']);

    fig1=figure('Visible','off'); imshow(BW_sauvola,[]); hold on;
    exportgraphics(fig1,binarisationPath,'BackgroundColor','none'); close(fig1);

    fig2=figure('Visible','off'); imshow(edgesCannyClean,[]); exportgraphics(fig2,edgePath,'BackgroundColor','none'); close(fig2);
    fig3=figure('Visible','off'); imshow(avgMonoCropped,[]); exportgraphics(fig3,avgPath,'BackgroundColor','none'); close(fig3);

    %% --- Plot flow detection ---
    fig4 = figure('Visible','off');
    plot(frameList, changeList, '-o', 'LineWidth', 1.2);
    hold on;

    % Baseline threshold
    bgWindow = 50; % same as in detectFlowFrames
    persistFrames = 3;
    clear yline
    if length(changeList) >= bgWindow
        baselineFrames = changeList(1:bgWindow);
        maxBg = max(baselineFrames);
        thresholdBaseline = maxBg * (1 + sigmaFactor);
        yline(thresholdBaseline, 'r--', 'LineWidth', 1.2);
        yline(maxBg, 'g--', 'LineWidth', 1.2);
    end

    xline(exactFlowFrame, 'b--', 'LineWidth', 1.5);

    xlabel('Frame Number');
    ylabel('Mean Frame-to-Frame Change');
    title(['Flow Detection for ', baseName]);
    legend('Mean Change','Baseline Threshold','Max of First Frames','Exact Start','Location','best');
    grid on;

    % Save figure
    flowPlotPath = fullfile(outputDir,[baseName '_FlowDetection.png']);
    exportgraphics(fig4, flowPlotPath, 'BackgroundColor','none');
    close(fig4);
    allResults(end+1,:) = {fileName, exactFlowFrame, sprayAngleLeft, sprayAngleRight, totalConeAngle, pixelSize_mm, ...
    binarisationPath, edgePath, avgPath, gifPath, flowPlotPath};

end

resultsTable = cell2table(allResults, ...
    'VariableNames',{'FileName','FirstFlowFrame','LeftSprayAngle_deg','RightSprayAngle_deg','TotalConeAngle_deg','PixelSize_mm', ...
                     'SauvolaMaskPath','CannyEdgesPath','AvgMonoImagePath','AllFramesGIFPath','FlowDetectionPlotPath'});

writetable(resultsTable,outputExcel);
fprintf('All results saved to Excel: %s\n', outputExcel);
