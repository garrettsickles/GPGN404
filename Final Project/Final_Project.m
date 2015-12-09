function Final_Project()
    % Get Filenames
    playedFile = uigetfile('*.*','Played audio file');
    recordFile = uigetfile('*.*','Recorded audio file');
    inputFile = uigetfile('*.*','Simulate audio file');
    
    % Import the data
    rawPlayed = importdata(playedFile);
    rawRecord = importdata(recordFile);
    rawInput = importdata(inputFile);
    
    % Check the data sampling frequencies
    disp('Sampling frequencies of the data')
    disp(['Played:   Fs = ', num2str(rawPlayed.fs)]);
    disp(['Recorded: Fs = ', num2str(rawRecord.fs)]);
    disp(['Simulate: Fs = ', num2str(rawInput.fs)]);
    
    % Fix the played audio
    i = 1;
    d = rawPlayed.data(i);
    while d == 0
        i = i + 1;
        d = rawPlayed.data(i);
    end
    rawPlayed.data = rawPlayed.data(i:end);
    
    % Fix the data sampling frequencies
    
    record = struct('data', rawRecord.data(:,1), 'fs', rawRecord.fs);
    
    [NP, DP] = rat(rawRecord.fs./rawPlayed.fs);
    played = struct('data', resample(rawPlayed.data(:,1), NP, DP), 'fs', rawRecord.fs);
    
    [NI, DI] = rat(rawRecord.fs./rawInput.fs);
    input = struct('data', resample(rawInput.data(:,1), NI, DI), 'fs', rawRecord.fs);
    
    % Deconvolve the signal
    RL = length(record.data);
    PL = length(played.data);
    SL = length(input.data);
    N = max(RL, PL, SL);
    thing1 = fft(record.data, N)./fft(played.data, N);
    filter = struct('data',ifft(fft(record.data, N)./fft(played.data, N)), 'fs', rawRecord.fs);
    
    % Plot the the filter
    figure('position', [0, 0, 800, 450]);
    subplot(2, 1, 1);
    plot(Tn(filter.data, filter.fs, 0.0), filter.data);
    title('Time Domain Response of the Filter');
    xlabel('Time (sec)');
    ylabel('Amplitude');
    axis tight;
    subplot(2, 1, 2);
    plot(Fk(Tn(filter.data, fliter.fs, 0.0)), thing1);
    title('Frequency Domain Response of the Filter');
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    axis tight;
    
    % Save the filter to a file
    audiowrite('filter.wav', filter.data, filter.fs);
    
    % Simulate a noise through the filter
    thing2 = fft(input.data, N).*fft(filter.data, N);
    output = struct('data', ifft(thing2), 'fs', filter.fs);
    
    % Plot the the sound
    figure('position', [0, 0, 800, 450]);
    subplot(2, 1, 1);
    plot(Tn(output.data, output.fs, 0.0), filter.data);
    title('Time Domain Response of the Simulated Noise');
    xlabel('Time (sec)');
    ylabel('Amplitude');
    axis tight;
    subplot(2, 1, 2);
    plot(Fk(Tn(output.data, output.fs, 0.0)), thing2);
    title('Frequency Domain Response of the Simulated Noise');
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    axis tight;
    
    % Write the simulated output
    audiowrite('Output.wav', output.data, output.fs);
end

% Produce a list of evenly spaced times given a dataset, Sampling Frequency, and initial time
function [result] = Tn(input, Fs, init)
    result = init + (0:1:(length(input)-1))*(1/Fs);
end

% Produce a list of frequencies given a list of evenly spaced times
function [result] = Fk(times)
    N = length(times);
    Ts = (times(N) - times(1)) / (N - 1);
    dF = 1 / (N * Ts);
    result = ((0:N-1) - ceil(N/2))*dF;
end