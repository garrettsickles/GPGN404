function Lab08()
    % Problem1();
    % Problem2();
    Problem3();
end

function Problem1()
    Fs = 25.0;
    Ts = 1/Fs;
    
    Tmin = 0.0;
    Tmax = 3.0;
    
    N = ceil((Tmax-Tmin)*Fs)
    n = 1:N;
    ts = (n-1).*Ts;
    tc = Tmin:0.001:Tmax;
    
    x = (2.0)*cos(3.1415926*(10*ts+1/6))+cos(3.1415926*(20*ts+1/3));  
    y = (2.0)*cos(3.1415926*(10*tc+1/6))+cos(3.1415926*(20*tc+1/3));
    z = zeros(length(tc),1);
    
    for i=1:length(tc)
        z(i) = sum(x .* sinc(tc(i)./Ts-(0:N-1)));
    end
    
    figure('position', [0, 0, 750, 250]);
    plot(tc, y, 'color', 'b');
    hold on;
    plot(tc, z, 'color', [1,0.3,0]);
    hold on;
    plot(ts, x, 'o', 'color', 'k');
    xlabel('Time (sec)');
    ylabel('Amplitude');
    legend('Original','Reconstructed','Samples Points');
    title('Problem 1: [0.0, 3.0] seconds');
end

function Problem2()
    Fs = 15.0;
    Ts = 1/Fs;
    
    Tmin = 0.0;
    Tmax = 3.0;
    
    N = ceil((Tmax-Tmin)*Fs)
    n = 1:N;
    ts = (n-1).*Ts;
    tc = Tmin:0.001:Tmax;
    
    x = (2.0)*cos(3.1415926*(10*ts+1/6))+cos(3.1415926*(20*ts+1/3));  
    y = (2.0)*cos(3.1415926*(10*tc+1/6))+cos(3.1415926*(20*tc+1/3));
    z = zeros(length(tc),1);
    
    for i=1:length(tc)
        z(i) = sum(x .* sinc(tc(i)./Ts-(0:N-1)));
    end
    
    figure('position', [0, 0, 750, 250]);
    plot(tc, y, 'color', 'b');
    hold on;
    plot(tc, z, 'color', [1,0.3,0]);
    hold on;
    plot(ts, x, 'o', 'color', 'k');
    xlabel('Time (sec)');
    ylabel('Amplitude');
    legend('Original','Reconstructed','Samples Points');
    title('Problem 2: [0.0, 3.0] seconds');
end

function [dates_interp, interpTmin, interpTmax] = Problem3()
    % Load the data
    filename = 'Lab6_t_T.csv';
    rawData = importdata(filename);
    rawData = rawData.data;
    rawDates = rawData(:, 1);
    % rawTmax = rawData(:, 2);
    rawTmin = rawData(:, 3);
    invalid = -9999;
    
    % Fill in the missing days
    serialDates = datenum(num2str(rawDates), 'yyyymmdd');
    serialDates = serialDates - serialDates(1);
    dates = (serialDates(1):1:serialDates(length(serialDates))) - serialDates(1);
    
    Tmin = fixData(serialDates, rawTmin, invalid);
    % Tmax = fixData(serialDates, rawTmax, invalid);
    
    % Interpolate the data
    Ts_interp = 4;
    dates_interp = dates(1):Ts_interp:dates(length(dates));
    interpTmin = sincInterp(dates, Tmin, dates_interp, Ts_interp, 1e-1);
    % interpTmax = sincInterp(dates, Tmax, dates_interp, Ts_interp, 1e-1);
    
    figure('position', [0, 0, 750, 450]);
    subplot(2,1,1);
    plot(dates, Tmin, 'color', 'b');
    axis tight;
    title('Minimum Temperature: Raw Data');
    xlabel('# of Days since 1 January 1893');
    ylabel('Tenths of ^{o}C');
    legend('Raw T_{min}');
    subplot(2,1,2);
    plot(dates_interp, interpTmin, 'color', [1,0.3,0]);
    axis tight;
    title('Minimum Temperature: Sinc Interpolated');
    xlabel('# of Days since 1 January 1893');
    ylabel('Tenths of ^{o}C');
    legend('Interpolated T_{min}');

%     figure('position', [0, 0, 750, 450]);
%     subplot(2,1,1);
%     plot(dates, Tmax, 'color', 'b');
%     axis tight;
%     title('Maximum Temperature: Raw Data');
%     xlabel('# of Days since 1 January 1893');
%     ylabel('Tenths of ^{o}C');
%     legend('Raw T_{max}');
%     subplot(2,1,2);
%     plot(dates_interp, interpTmax, 'color', [1,0.3,0]);
%     axis tight;
%     title('Maximum Temperature: Sinc Interpolated');
%     xlabel('# of Days since 1 January 1893');
%     ylabel('Tenths of ^{o}C');
%     legend('Interpolated T_{max}');
    
    Tmin_Xk = fftshift(fft(Tmin))/length(Tmin);
%    Tmax_Xk = fftshift(fft(Tmax))/length(Tmax);
    Fk = getFrequencies(dates.*(24.0*60.0*60.0));
    
    Ff = 1/(365.25*24.0*60.0*60.0);
    fTmin_Xk(abs(Tmin_Xk) < Ff) = 0.0;
%    fTmax_Xk(abs(Tmax_Xk) < Ff) = 0.0;
    
    invTmin = ifft(ifftshift(fTmin_Xk));
%    invTmax = ifft(ifftshift(fTmax_Xk));
    
    figure('position', [0, 0, 750, 450]);
    subplot(2,1,1);
    plot(Fk, abs(Tmin_Xk));
    axis tight;
    xlabel('Frequncy (Hz)');
    ylabel('Amplitude');
    title('Amplitude Spectrum of Raw T_{min}');
    legend('Raw T_{min}');
    subplot(2,1,2);
    plot(Fk, abs(fTmin_Xk));
    axis tight;
    xlabel('Frequncy (Hz)');
    ylabel('Amplitude');
    title('Amplitude Spectrum of Filtered T_{min}');
    legend('Filtered T_{min}');
    
%     figure('position', [0, 0, 750, 450]);
%     subplot(2,1,1);
%     plot(Fk, abs(Tmax_Xk));
%     axis tight;
%     xlabel('Frequncy (Hz)');
%     ylabel('Amplitude');
%     title('Amplitude Spectrum of Raw T_{max}');
%     legend('Raw T_{max}');
%     subplot(2,1,2);
%     plot(Fk, abs(fTmax_Xk));
%     axis tight;
%     xlabel('Frequncy (Hz)');
%     ylabel('Amplitude');
%     title('Amplitude Spectrum of Filtered T_{max}');
%     legend('Filtered T_{max}');
    
    figure('position', [0, 0, 750, 450]);
    subplot(2,1,1);
    plot(dates, Tmin, 'color', 'b');
    axis tight;
    title('Minimum Temperature: Raw Data');
    xlabel('# of Days since 1 January 1893');
    ylabel('Tenths of ^{o}C');
    legend('Raw T_{min}');
    subplot(2,1,2);
    plot(dates_interp, invTmin, 'color', [1,0.3,0]);
    axis tight;
    title('Minimum Temperature: Interpolated Filtered Data');
    xlabel('# of Days since 1 January 1893');
    ylabel('Tenths of ^{o}C');
    legend('Filtered T_{min}');
    
%     figure('position', [0, 0, 750, 450]);
%     subplot(2,1,1);
%     plot(dates, Tmin, 'color', 'b');
%     axis tight;
%     title('Maximum Temperature: Raw Data');
%     xlabel('# of Days since 1 January 1893');
%     ylabel('Tenths of ^{o}C');
%     legend('Raw T_{max}');
%     subplot(2,1,2);
%     plot(dates_interp, invTmin, 'color', [1,0.3,0]);
%     axis tight;
%     title('Maximum Temperature: Interpolated Filtered Data');
%     xlabel('# of Days since 1 January 1893');
%     ylabel('Tenths of ^{o}C');
%     legend('Filtered T_{max}');
end

% Fill in missing data
function [result] = fixData(input, values, invalid)
    result(1) = values(1);
    for i=2:length(input)
        d = input(i) - input(i-1);
        if d > 1
            vertcat(result, nan(d-1, 1));
        end
        vertcat(result, values(i));
    end
    result(isnan(result)) = [];
    result(result == invalid) = [];
end

% Sinc Interpolate the data, by Jyoti Behura
function [result] = sincInterp(input, values, output, Ts, tolerance)
    ni = length(input);
    no = length(output);
    A = zeros(ni,no);
    for i=0:ni-1
        for j=0:no-1
            A(i+1,j+1) = sinc((input(i+1)-output(j+1))/Ts);
        end
    end
    result = pinv(A, tolerance)*values;
end

% Produce a list of frequencies given a list of evenly spaced times
function [result] = getFrequencies(times)
    N = length(times);
    Ts = (times(N) - times(1)) / (N - 1);
    dF = 1 / (N * Ts);
    result = ((0:N-1) - ceil(N/2))*dF;
end