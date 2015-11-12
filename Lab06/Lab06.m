% --------------------- Question 1 -----------------------
function Lab06()
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
    
    for i=2:length(fixedDates)
        if fixedTmax(i) == invalid
            fixedTmax(i) = fixedTmax(i-1);
        end
        if fixedTmin(i) == invalid
            fixedTmin(i) = fixedTmin(i-1);
        end
    end
    
    
    MAF = MovingAverage(fixedTmax, period);
    
    figure;
    % plot(serialDates, [Tmax, MAF((period-1)/2:(period-1)/2 - 1 + length(serialDates))]);
    plot(fixedDates, fixedTmax);
    hold on;
    plot(MAF(1:length(MAF)));
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
    result = MovingAverage(input);
    for i=1:number
        result(i) = result(i) * (number / i);
    end
end
