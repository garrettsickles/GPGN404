% --------------------- Question 1 -----------------------
function Lab06()
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
    
    for i=2:length(dates)
        if Tmax(i) == -9999
            Tmax(i) = Tmax(i-1);
        end
        if Tmin(i) == -9999
            Tmin(i) = Tmin(i-1);
        end
    end
    MAF = MovingAverage(Tmax, period);
    
    figure;
    plot(serialDates, [Tmax, MAF((period-1)/2:(period-1)/2 - 1 + length(serialDates))]);
    
    figure;
    plot(serialDates);
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
end
