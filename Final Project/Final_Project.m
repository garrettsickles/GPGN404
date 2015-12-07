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
    disp('Number of data points')
    disp(['Played:   ', num2str(length(record.data))]);
    disp(['Recorded: ', num2str(length(played.data))]);
    disp(['Simulate: ', num2str(length(input.data))]);
    fdata = deconv(record.data, played.data);
    filter = struct('data', fdata, 'fs', rawRecord.fs);
    audiowrite('filter.wav', filter.data, filter.fs);
    output = struct('data', conv(input.data, filter.data), 'fs', filter.fs);
    
    % Write the simulated output
    audiowrite('Output.wav', output.data, output.fs);
end