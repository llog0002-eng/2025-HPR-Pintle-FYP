function exactFlowFrame = detectFlowFrames(filePath, monoChannel, startFrame, maxFrames, sigmaFactor, maxPixelValue)
    % DETECTFLOWFRAMES - Detect first flow frame using robust background estimation
    % Inputs:
    % filePath      : path to raw video
    % monoChannel   : 1=R, 2=G, 3=B
    % startFrame    : frame to start searching
    % maxFrames     : maximum frames to check
    % sigmaFactor   : multiplier for robust std (increase to reduce sensitivity)
    % maxPixelValue : maximum pixel value (e.g., 4095 for 12-bit)

    % --- Preallocate storage ---
    frameList = [];
    changeList = [];

    % --- Read the first frame ---
    prevFrameMono = double(readmraw(filePath, startFrame));
    prevFrameMono = prevFrameMono(:,:,monoChannel);

    % --- Compute frame-to-frame mean changes ---
    for f = startFrame+1:maxFrames
        I = double(readmraw(filePath, f));
        I_mono = I(:,:,monoChannel);

        meanChange = mean(abs(I_mono - prevFrameMono), 'all') / maxPixelValue;

        frameList(end+1) = f;
        changeList(end+1) = meanChange;

        prevFrameMono = I_mono;
    end

    % --- Robust background estimation ---
    medianBg = median(changeList);                  % Median of all frame changes
    madBg = median(abs(changeList - medianBg));    % Median absolute deviation
    thresholdDynamic = medianBg + sigmaFactor * madBg;

    % --- Detect first frame exceeding threshold ---
    idx = find(changeList > thresholdDynamic, 1, 'first');
    if isempty(idx)
        exactFlowFrame = startFrame;
        disp('No flow detected.');
    else
        exactFlowFrame = frameList(idx);
        disp(['Exact first flow frame detected: ', num2str(exactFlowFrame)]);
    end

    % --- Plot results ---
    figure;
    plot(frameList, changeList, '-o');
    hold on;
    yline(thresholdDynamic, 'r--', 'Dynamic Threshold');
    xline(exactFlowFrame, 'b--', 'Exact Flow Start');
    xlabel('Frame Number');
    ylabel('Mean Frame-to-Frame Change');
    title('Flow Detection over Time');
    legend('Mean Change','Dynamic Threshold','Exact Start');
    grid on;
end
