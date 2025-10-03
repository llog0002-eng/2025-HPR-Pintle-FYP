function exactFlowFrame = detectFlowFrames(filePath, monoChannel, startFrame, maxFrames, sigmaFactor, maxPixelValue)
% DETECTFLOWFRAMES - Detect first flow frame using max of initial frames as threshold
% Automatically reads more frames if flow not detected in initial maxFrames
%
% Inputs:
%   filePath      : path to raw video
%   monoChannel   : 1=R, 2=G, 3=B
%   startFrame    : frame to start searching
%   maxFrames     : initial number of frames to check
%   sigmaFactor   : multiplier to slightly raise threshold above max
%   maxPixelValue : maximum pixel value (e.g., 4095 for 12-bit)

    % --- Preallocate ---
    frameList = [];
    changeList = [];

    % --- Read first frame ---
    prevFrameMono = double(readmraw(filePath, startFrame));
    prevFrameMono = prevFrameMono(:,:,monoChannel);

    % --- Parameters ---
    bgWindow = 50;        % Number of initial frames to estimate baseline
    persistFrames = 3;    % Require this many consecutive frames above threshold

    exactFlowFrame = startFrame;
    flowDetected = false;

    f = startFrame + 1;   % Current frame index

    while ~flowDetected
        % --- Read frame ---
        I = double(readmraw(filePath, f));
        I_mono = I(:,:,monoChannel);

        % --- Compute mean frame-to-frame change ---
        meanChange = mean(abs(I_mono - prevFrameMono), 'all') / maxPixelValue;

        frameList(end+1) = f;
        changeList(end+1) = meanChange;

        prevFrameMono = I_mono;

        % --- Baseline threshold from first bgWindow frames ---
        if length(changeList) >= bgWindow && ~exist('thresholdBaseline', 'var')
            baselineFrames = changeList(1:bgWindow);
            maxBg = max(baselineFrames);
            thresholdBaseline = maxBg * (1 + sigmaFactor);  % Slightly above max to allow drift
        end

        % --- Detect flow ---
        if exist('thresholdBaseline', 'var') && length(changeList) >= persistFrames
            recentChange = changeList(end-persistFrames+1:end);
            if all(recentChange > thresholdBaseline)
                exactFlowFrame = f - persistFrames + 1;
                flowDetected = true;
                break;
            end
        end

        f = f + 1;

        % --- Optionally, you could check if f exceeds total frames in video ---
        % If your video has a fixed length, stop when end is reached
        % For example, totalFrames = getTotalFrames(filePath);
        % if f > totalFrames, break
    end

    if ~flowDetected
        disp('No flow detected within frames read so far.');
    else
        disp(['Exact first flow frame detected: ', num2str(exactFlowFrame)]);
    end

    % --- Plot ---
    figure;
    plot(frameList, changeList, '-o');
    hold on;

    if exist('thresholdBaseline', 'var')
        yline(thresholdBaseline, 'r--', 'Baseline Threshold');
        yline(maxBg, 'g--', 'Max of First Frames');  % Show raw max
    end
    
    xline(exactFlowFrame, 'b--', 'Exact Flow Start');

    xlabel('Frame Number');
    ylabel('Mean Frame-to-Frame Change');
    title('Flow Detection with Max Initial Baseline Threshold');
    legend('Mean Change','Baseline Threshold','Max of First Frames','Exact Start');
    grid on;
end
