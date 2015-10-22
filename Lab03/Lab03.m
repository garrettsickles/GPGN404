%------------- Problem 1 -----------------
Fk = [0 5 20];
Xk = [5 1.5*exp(1i*pi*0.25) 2.0*exp(1i*pi*0.5)];
Amp = Spectra(Fk, Xk, 'Amplitude');
Phs = Spectra(Fk, Xk, 'Phase');
clear;

%------------- Problem 3 -----------------

figure;
t3 = -1:0.001:1;
plot(t3, 5 + 3.0*cos(2*3.14159*5.0*t3+3.14159/4.0)+4.0*cos(2*3.14159*20.0*t3+3.14159/2.0));
xlabel({'Time (sec)'});
title('Time Continuous Signal');
ylabel({'Amplitude'});
clear;

%------------- Problem 4 -----------------
Fs = 100;            % Sampling frequency
T = 1/Fs;             % Sampling period
L = 1000;             % Length of signal
t = (0:L-1)*T;        % Time vector
X = 5 + 3.0*cos(2*3.14159*5.0*t+3.14159/4.0)+4.0*cos(2*3.14159*20.0*t+3.14159/2.0);
Y = fftshift(fft(X))/(length(t));
f = Fs*(-L/2:L/2-1)/L;
Phase = zeros(length(Y), 1);
Amplitude = abs(Y);
for i = 1:L
    Phase(i) = angle(Y(i));
    %
    % if Amplitude(i) <= 0.001
    %    Amplitude(i) = NaN;
    %    Phase(i) = NaN;
    % end
end
figure;
stem(f,Amplitude);
xlabel({'Frequency (Hz)'});
title('FFT Amplitude Spectrum');
ylabel({'Amplitude'});
figure;
stem(f,Phase);
xlabel({'Frequency (Hz)'});
title('FFT Phase Spectrum');
ylabel({'Phase'});
clear;

%------------- Problem 5 -----------------
figure;
Data = load('Lab3_t_xt.dat');
TimeData = Data(:,1);
AmpData = Data(:,2);
plot(TimeData, AmpData);
xlabel({'Time (sec)'});
title('Time Continuous Signal from Data file');
ylabel({'Amplitude'});

Fs = 5000;           % Sampling frequency
T = 1/Fs;             % Sampling period
t = TimeData;       % Time vector

X = AmpData;
Y = fftshift(fft(X))/(length(t));
fk = Fs/2*linspace(-1,1,length(t));
Phase = angle(Y);
Amplitude = abs(Y);

figure;
stem(fk,Amplitude);
xlabel({'Frequency (Hz)'});
title('Data Amplitude Spectrum');
ylabel({'Amplitude'});
figure;
stem(fk,Phase);
xlabel({'Frequency (Hz)'});
title('Data Phase Spectrum');
ylabel({'Phase'});
clear;