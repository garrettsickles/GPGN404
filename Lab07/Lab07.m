function Lab07()
    filename = uigetfile('*.wav','Input WAV audio file');
    rawData = importdata(filename);
    Fs = rawData.fs;
    data = rawData.data;
    init = 0.0;
    fps = 10;
    
    times = getTimes(data, Fs, 0.0);
    frameCount = ceil((times(length(times))-times(1))*fps);
    pointsPerFrame = ceil(length(times)/frameCount);
    min = 0.0;
    max = 0.25;
    
    windowedTimes = getWindowedTimes(data, Fs, fps, init);
    windowedData = getWindowedData(data, Fs, fps, init);
    
    FFTwindowedData = fftshift(fft(windowedData,pointsPerFrame, 2),2).*(1/pointsPerFrame);
    fk = getFrequencies(windowedTimes(1,:));
    
    times = PadTimes(times, frameCount*pointsPerFrame);
    data = PadData(data, frameCount*pointsPerFrame);
    v = VideoWriter('animation.avi');
    
    v.FrameRate = fps;
    open(v);
    figure;
    RGBblue = [51 102 255]/255;
    RGBgreen = [0 255 0]/255;
    for i=1:size(windowedData, 1)
        subplot(5,1,1:3);
        plot(fk, abs(FFTwindowedData(i,:)), 'Color', RGBblue);
        title(['Spectrogram of Audio File']);
        xlabel('Frequency (Hz)');
        ylabel('Amplitude');
        axis tight;
        ylim([min max]);
        subplot(5,1,5);
        plot(times, data, 'Color', RGBblue);
        hold on;
        plot(windowedTimes(i,:), windowedData(i,:), 'Color', RGBgreen);
        title(['Windowed Time Domain Data'])
        xlabel('Time (sec)');
        ylabel('Amplitude');
        axis tight;
        writeVideo(v, getframe(gcf));
    end
    close(v);
end

function [result] = getWindowedData(input, Fs, fps, init)
    times = getTimes(input, Fs, init);
    frameCount = ceil((times(length(times))-times(1))*fps);
    pointsPerFrame = ceil(length(input)/frameCount);
    padded = PadData(input, frameCount*pointsPerFrame);
    result = reshape(padded, pointsPerFrame, frameCount)';
end

function [result] = getWindowedTimes(input, Fs, fps, init)
    times = getTimes(input, Fs, init);
    frameCount = ceil((times(length(times))-times(1))*fps);
    pointsPerFrame = ceil(length(input)/frameCount);
    times = PadTimes(times, frameCount*pointsPerFrame);
    result = reshape(times, pointsPerFrame, frameCount)';
end

% Pad the data using a reverse cosine taper to maintain frequencies
function [result] = PadData(input, l)
    lengthInput = length(input);
    padLength = l - lengthInput;
    result = vertcat(input, zeros(padLength, 1));
    if padLength > 1
        dx = (3.1415926/2.0)/(padLength-1);
        x = 0.0:dx:(3.1415926/2.0);
        reverse = flipud(result(lengthInput-padLength+1:lengthInput));
        reverse = ((reverse - input(lengthInput)).*-1)+input(lengthInput);
        result((lengthInput+1):length(result)) = reverse'.*cos(x);
    end
end

% Pad the times for a matrix
function [result] = PadTimes(times, l)
    padLength = l - length(times);
    result = times;
    if padLength > 0
        dt = times(2)-times(1);
        addTimes = (1:1:padLength).*dt+times(length(times));
        result = horzcat(result, addTimes);
    end
end

% Produce a list of evenly spaced times given a dataset, Sampling Frequency, and initial time
function [result] = getTimes(input, Fs, init)
    result = init + (0:1:(length(input)-1))*(1/Fs);
end

% Produce a list of frequencies given a list of evenly spaced times
function [result] = getFrequencies(times)
    N = length(times);
    Ts = (times(N) - times(1)) / (N - 1);
    dF = 1 / (N * Ts);
    result = ((0:N-1) - ceil(N/2))*dF;
end