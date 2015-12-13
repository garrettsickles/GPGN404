function Final_Project(title)
    % Import the data
    junk1 = importdata(uigetfile('*.*','Played audio file'));
    junk2 = importdata(uigetfile('*.*','Recorded audio file'));
    junk3 = importdata(uigetfile('*.*','Simulate audio file'));
    rawPlayed = junk1.data(:,1);
    RPFS = junk1.fs;
    rawRecord = junk2.data(:,1);
    RRFS = junk2.fs;
    rawInput = junk3.data(:,1);
    RIFS = junk3.fs;
    % Fix the played audio
    i = 1;
    d = rawPlayed(i);
    while d == 0
        i = i + 1;
        d = rawPlayed(i);
    end
    rawPlayed = rawPlayed(i:end);
    % Fix the data sampling frequencies
    Fs = RPFS;
    [NP, DP] = rat(Fs./RPFS);
    played = struct('data', resample(rawPlayed(:,1), NP, DP), 'fs', Fs);
    [NI, DI] = rat(Fs./RIFS);
    input = struct('data', resample(rawInput(:,1), NI, DI), 'fs', Fs);
    [NR, DR] = rat(Fs./RRFS);
    record = struct('data', resample(rawRecord(:,1), NR, DR), 'fs', Fs); 
    % Deconvolve the signal
    N = max([length(record.data) length(played.data) length(input.data)]);
    filter = struct('data',ifft(fft(record.data, N)./fft(played.data, N)), 'fs', record.fs);
    % Band Pass 
    filter.data = bandPass(filter.data, filter.fs, 20.0, 2000.0);
    % Plot the the played
    time2freq(played.data, played.fs, 'Played Frequency Sweep');
    % Plot the the recorder
    time2freq(record.data, record.fs, ['Recorded in ', title]);
    % Plot the the input
    time2freq(input.data, input.fs, 'Input: Mr. Mackey');
    % Plot the the filter
    time2freq(filter.data, filter.fs, ['Filter in ', title]);
    % Simulate a noise through the filter
    output = struct('data', ifft(fft(input.data, N).*fft(filter.data, N)), 'fs', filter.fs);
    % Plot the the sound
    output.data = output.data ./ max(output.data);
    time2freq(output.data, output.fs, [' Mr. Mackey in ', title]);
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

function [result] = bandPass(input, fs, fmin, fmax)
    freq = Fk(Tn(input, fs, 0.0));
    filterFFT = fftshift(fft(input));
    filterFFT(abs(freq) < fmin) = 0.0;
    filterFFT(abs(freq) > fmax) = 0.0;
    result = ifft(fftshift(filterFFT));
end

% Plot a time and frequency response of a data set
function [result] = time2freq(input, fs, name)
    result = fftshift(fft(input));
    time = Tn(input, fs, 0.0);
    figure('position', [0, 0, 800, 450]);
    subplot(2, 1, 1);
    plot(time, input);
    title(['Time Domain Response of the ', name]);
    xlabel('Time (sec)');
    ylabel('Amplitude');
    axis tight;
    subplot(2, 1, 2);
    plot(Fk(time), result);
    title(['Frequency Domain Response of the ', name]);
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    axis tight;
end