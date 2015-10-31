function Lab05()
    clear;

    filename = 'data.dat';
    data = tblread(filename);
    time = data(:,1);
    value = data(:,2);
    trans = DFT(value, 'forward');
    inverse = DFT(trans, 'inverse');
    
    figure('position', [0, 0, 700, 600]);
    orient tall;
    subplot(5,1,1);
    plot(time, value);
    axis tight;
    xlabel('Time (sec)');
    ylabel('Amplitude');
    title('Time Signal from file');
    
    subplot(5,1,2);
    stem(Fk(time), abs(trans));
    axis tight;
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    title('Amplitude Spectrum using D.F.T.');
    
    subplot(5,1,3);
    stem(Fk(time), angle(trans));
    axis tight;
    xlabel('Frequency (Hz)');
    ylabel('Phase (rad)');
    title('Phase Spectrum using D.F.T.');
    
    subplot(5,1,4);
    plot(time, inverse);
    axis tight;
    xlabel('Time (sec)');
    ylabel('Amplitude');
    title('Calculated signal using inverse D.F.T.');
    
    subplot(5,1,5);
    plot(time, value - inverse);
    axis tight;
    xlabel('Time (sec)');
    ylabel('Amplitude');
    title('Calculated difference of inverted and actual signal');
end

function [result] = Fk(times)
    N = length(times);
    Ts = (times(N) - times(1)) / (N - 1);
    dF = 1 / (N * Ts);
    result = ((0:N-1) - ceil(N/2))*dF;
end

function [result] = DFT(values, direction)
    N = length(values);
    d = (strcmpi(direction, 'forward') * 2) - 1;
    matrix = zeros(length(values));
    Nyquist = ceil(N / 2);
    if d == -1
        values = vertcat(values(Nyquist+1:N), values(1:Nyquist));
    end
    for n = 0:N-1
        for k = 0:n
            if n == 0
                operand = 1;
            else
                operand = exp(2i*3.1415926*n*k*d/N);
            end
            matrix(n+1,k+1) = operand;
            if n ~= k
                matrix(k+1,n+1) = operand;
            end
        end
    end
    result = matrix * values;
    if d == 1
        result = result .* (1/N);
        result = vertcat(result(Nyquist+1:N),result(1:Nyquist));
    end
end