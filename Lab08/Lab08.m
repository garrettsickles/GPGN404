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

function Problem3()
    invalid = -9999;
    filename = 'Lab6_t_T.csv';
    data = importdata(filename);
    data = data.data;
    dates = data(:, 1);
    Tmax = data(:, 2);
    Tmin = data(:, 3);
    period = 365;
    
    serialDates = datenum(num2str(dates), 'yyyymmdd');
    serialDates = serialDates - serialDates(1);
    
    fixedDates = serialDates(1):1:serialDates(length(serialDates));
    N = length(fixedDates);
    
    fixedTmin(1) = Tmin(1);
    fixedTmax(1) = Tmax(1);
    for i=2:length(serialDates)
        d = serialDates(i) - serialDates(i-1);
        if d > 1
            fixedTmax = vertcat(fixedTmax, ones(d-1, 1).*invalid); 
            fixedTmin = vertcat(fixedTmin, ones(d-1, 1).*invalid);
        end
        fixedTmax = vertcat(fixedTmax, Tmax(i));
        fixedTmin = vertcat(fixedTmin, Tmin(i));
    end
    for i=1:length(fixedDates)
        if fixedTmax(i) == invalid
            fixedTmax(i) = NaN;
        end
        if fixedTmin(i) == invalid
            fixedTmin(i) = NaN;
        end
    end
    rawTmin = fixedTmin;
    rawTmax = fixedTmax;
    for i=1:length(fixedDates)*4
        %if isnan(fixedTmax(i))
            fixedTmax(i) = nansum(rawTmax' .* sinc(fixedDates(i)/4-(0:N-1)));
        %end
        %if isnan(fixedTmin(i))
            fixedTmin(i) = nansum(rawTmin' .* sinc(fixedDates(i)/4-(0:N-1)));
        %end
    end
    
    
    MAF_Tmax = MovingAverage(fixedTmax, period);
    NLMAF_Tmax = NoLagMovingAverage(fixedTmax, period);
    MAF_Tmin = MovingAverage(fixedTmin, period);
    NLMAF_Tmin = NoLagMovingAverage(fixedTmin, period);
    
    figure('position', [0, 0, 750, 450]);
    subplot(2,1,1);
    plot(fixedDates, rawTmin, 'color', 'b');
    axis tight;
    title('Minimum Temperature: Raw Data');
    xlabel('# of Days since 1 January 1893');
    ylabel('Tenths of ^{o}C');
    legend('Raw T_{min}');
    subplot(2,1,2);
    plot(fixedDates, fixedTmin, 'color', [1,0.3,0]);
    axis tight;
    title('Minimum Temperature: Sinc Interpolated');
    xlabel('# of Days since 1 January 1893');
    ylabel('Tenths of ^{o}C');
    legend('Interpolated T_{min}');

    figure('position', [0, 0, 750, 450]);
    subplot(2,1,1);
    plot(fixedDates, rawTmax, 'color', 'b');
    axis tight;
    title('Maximum Temperature: Raw Data');
    xlabel('# of Days since 1 January 1893');
    ylabel('Tenths of ^{o}C');
    legend('Raw T_{max}');
    subplot(2,1,2);
    plot(fixedDates, fixedTmax, 'color', [1,0.3,0]);
    axis tight;
    title('Maximum Temperature: Sinc Interpolated');
    xlabel('# of Days since 1 January 1893');
    ylabel('Tenths of ^{o}C');
    legend('Interpolated T_{max}');
    
%     figure('position', [0, 0, 600, 250])
%     
%     figure('position', [0, 0, 600, 250]);
%     plot(fixedDates, fixedTmax, 'color', [0.0,0.0,0.0]+0.6);
%     hold on;
%     plot(NLMAF_Tmax(1:length(NLMAF_Tmax)));
%     axis tight;
%     title('Maximum Temperature No Lag Moving Average');
%     xlabel('# of Days since 1 January 1893');
%     ylabel('Tenths of ^{o}C');
%     legend('T_{max}','T_{max} Moving Average');
%     
%     figure('position', [0, 0, 600, 250]);
%     plot(fixedDates, fixedTmin, 'color', [0.0,0.0,0.0]+0.6);
%     hold on;
%     plot(NLMAF_Tmin(1:length(NLMAF_Tmin)));
%     axis tight;
%     title('Minimum Temperature No Lag Moving Average');
%     xlabel('# of Days since 1 January 1893');
%     ylabel('Tenths of ^{o}C');
%     legend('T_{min}','T_{min} Moving Average');
%     
%     figure('position', [0, 0, 600, 250]);
%     subplot(1,2,1);
%     plot(fixedDates, fixedTmax, 'color', [0.0,0.0,0.0]+0.6);
%     hold on;
%     plot(MAF_Tmax(1:length(MAF_Tmax)));
%     axis tight;
%     title('Maximum Temperature Raw Moving Average');
%     xlabel('# of Days since 1 January 1893');
%     ylabel('Tenths of ^{o}C');
%     legend('T_{max}','T_{max} Moving Average');
%     
%     subplot(1,2,2);
%     plot(fixedDates, fixedTmax, 'color', [0.0,0.0,0.0]+0.6);
%     hold on;
%     plot(NLMAF_Tmax(1:length(NLMAF_Tmax)));
%     axis tight;
%     title('Maximum Temperature No Lag Moving Average');
%     xlabel('# of Days since 1 January 1893');
%     ylabel('Tenths of ^{o}C');
%     legend('T_{max}','T_{max} Moving Average');
end

% Discrete Convolution 
function [result] = convolve(input, impulse)
    result = zeros(length(input) + length(impulse) - 1,1);
    for k=1:length(input)
        for j=1:length(impulse)
            result(k + j - 1) = result(k + j - 1) + impulse(j) * input(k);
        end
    end
end

% Moving Average Filter
function [result] = MovingAverage(input, number)
    impulse = ones(number, 1) ./ number;
    result = convolve(input, impulse);
end

% Tapered Moving Average Filter
function [result] = NoLagMovingAverage(input, number)
    result = MovingAverage(input, number);
    for i=1:number
        result(i) = result(i) * (number / i);
    end
end
