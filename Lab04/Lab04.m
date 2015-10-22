clear;

filename = 'data.dat';
data = tblread(filename);
time = data(:,1);
value = data(:,2);

% ----------------- 1a -----------------
%Fs = length(time) / time(length(time));
Fs = 1 / 0.020202;

Xk = fftshift(fft(value)) / length(time);
Fk = -Fs/2:Fs/length(time):Fs/2;
Fk = Fk(1:length(time));

figure('position', [0, 0, 750, 250]);
plot(time, value);
axis tight;
xlabel('Time (sec)');
ylabel('Amplitude');
title(['Time Signal: ',filename]);

figure('position', [0, 0, 750, 250]);
stem(Fk, abs(Xk));
axis tight;
xlabel('Frequncy (Hz)');
ylabel('Amplitude');
title(['Amplitude Spectrum: ',filename]);

DC = 0.0;
for i = 1:length(Fk)
    if Fk(i) == 0.0 || (Fk(i+1) > 0 && Fk(i) < 0)
        DC = Xk(i+1);
        break;
    end
end

% ----------------- 2b -----------------
Fs = 20;
Ts = 1.0 / Fs;

t1 = 0.0;
t2 = 2.0;
time = t1:Ts:t2;
if time(length(time)) >= 2.0
    time = time(1:(length(time)-1));
end

Xt = sin(10*3.1415926*time);

figure('position', [0, 0, 750, 250]);
stem(time, Xt);
axis tight;
xlabel('Time (sec)');
ylabel('Amplitude');
title(['Discrete Time Signal: sin(10\pit) sampled at ', num2str(Fs), ' hz']);

Xk = fftshift(fft(Xt)) / length(time);
Fk = -Fs/2:Fs/length(time):Fs/2;
Fk = Fk(1:length(time));

figure('position', [0, 0, 750, 250]);
stem(Fk, abs(Xk));
axis tight;
xlabel('Frequncy (Hz)');
ylabel('Amplitude');
title('Amplitude Spectrum: #2b');

time = -2.0:0.01:2.0;
Xr = zeros(1,length(time));
Xt = sin(10*3.1415926*time);
for i=length(Fk)/2:length(Fk)
    Xr = Xr + Xk(i) * exp(2*1i*3.14159*Fk(i)*time);
    Xr = Xr + conj(Xk(i)) * exp(-2*1i*3.14159*Fk(i)*time);
end

figure('position', [0, 0, 750, 750]);
subplot(3,1,1);
plot(time, Xr);
axis tight;
xlabel('Time (sec)');
ylabel('Amplitude');
title('Reconstructed Signal: #2b');
subplot(3,1,2);
plot(time, Xt);
axis tight;
xlabel('Time (sec)');
ylabel('Amplitude');
title('Actual Signal: #2b');
subplot(3,1,3);
plot(time, Xt-Xr);
axis tight;
xlabel('Time (sec)');
ylabel('Amplitude');
title('Signal Difference: #2b');

% ----------------- 2c -----------------
Fs = 9;
Ts = 1.0 / Fs;

t1 = 0.0;
t2 = 2.0;
time = t1:Ts:t2;
if time(length(time)) >= 2.0
    time = time(1:(length(time)-1));
end

Xt = sin(10*3.1415926*time);

figure('position', [0, 0, 750, 250]);
stem(time, Xt);
axis tight;
xlabel('Time (sec)');
ylabel('Amplitude');
title(['Discrete Time Signal: sin(10\pit) sampled at ', num2str(Fs), ' hz']);

Xk = fftshift(fft(Xt)) / length(time);
Fk = -Fs/2:Fs/length(time):Fs/2;
Fk = Fk(1:length(time));

figure('position', [0, 0, 750, 250]);
stem(Fk, abs(Xk));
axis tight;
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Amplitude Spectrum #2c');

time = -2.0:0.01:2.0;
Xr = zeros(1,length(time));
Xt = sin(10*3.1415926*time);

for i=length(Fk)/2:length(Fk)
    Xr = Xr + Xk(i) * exp(2*1i*3.14159*Fk(i)*time);
    Xr = Xr + conj(Xk(i)) * exp(-2*1i*3.14159*Fk(i)*time);
end

figure('position', [0, 0, 750, 750]);
subplot(3,1,1);
plot(time, Xr);
axis tight;
xlabel('Time (sec)');
ylabel('Amplitude');
title('Reconstructed Signal: #2c');
subplot(3,1,2);
plot(time, Xt);
axis tight;
xlabel('Time (sec)');
ylabel('Amplitude');
title('Actual Signal: #2c');
subplot(3,1,3);
plot(time, Xt-Xr);
axis tight;
xlabel('Time (sec)');
ylabel('Amplitude');
title('Signal Difference: #2c');

% ----------------- 3b -----------------
Fs = 100;
period = 1.0;
time = 0.0:1/Fs:4.0;

figure('position', [0, 0, 750, 250]);
plot(time, mod(time, period) < period*0.5);
axis tight;
xlabel('Time (sec)');
ylabel('Amplitude');
title(['Boxcar Signal Sampled at ', num2str(Fs), ' Hz: #3b']);

% ----------------- 3c + 3d -----------------
harmonic = 20;
Fs = 100;
period = 1.0;
time = 0.0:1/Fs:4.0;

Xt = zeros(1,length(time));

for k=-harmonic:harmonic
    if mod(k,2)
        Xt = Xt + ((1i*3.14159*k).^(-1)*exp((2i*3.14159*k*time)/period));
    end
end

Xt = Xt + 0.5;

Xk = fftshift(fft(Xt)) / length(time);
Fk = -Fs/2:Fs/length(time):Fs/2;
Fk = Fk(1:length(time));

figure('position', [0, 0, 750, 250]);
stem(Fk, abs(Xk));
axis tight;
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Amplitude Spectrum: #3c');

figure('position', [0, 0, 750, 250]);
stem(Fk, angle(Xk));
axis tight;
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');
title('Phase Spectrum: #3d');

% ----------------- 3e -----------------
harmonic = 100;
Fs = 100;
period = 1.0;
time = 0.0:1/Fs:4.0;

harmonic_matrix = zeros(2*harmonic+1, length(time));
for k=-harmonic:harmonic
    row = k + 1 + harmonic;
    if mod(k,2)
       harmonic_matrix(row,:) = (1i*3.14159*k).^(-1)*exp((2i*3.14159*k*time)/period);
    elseif k == 0
       harmonic_matrix(row,:) = harmonic_matrix(row) + 0.5;
    else
       harmonic_matrix(row,:) = 0.0;
    end
end

h = 3;
Xt = zeros(1,length(time));
for k=-h:h
    row = k + 1 + harmonic;
    Xt = Xt + harmonic_matrix(row,:);
end
figure('position', [0, 0, 750, 250]);
plot(time,Xt);
axis tight;
xlabel('Time (sec)');
ylabel('Amplitude');
title(['Synthesized Square Wave for ', num2str(h),' Harmonics Sampled at ', num2str(Fs), ' Hz: #3e']);

h = 10;
Xt = zeros(1,length(time));
for k=-h:h
    row = k + 1 + harmonic;
    Xt = Xt + harmonic_matrix(row,:);
end
figure('position', [0, 0, 750, 250]);
plot(time,Xt);
axis tight;
xlabel('Time (sec)');
ylabel('Amplitude');
title(['Synthesized Square Wave for ', num2str(h),' Harmonics Sampled at ', num2str(Fs), ' Hz: #3f']);

h = 20;
Xt = zeros(1,length(time));
for k=-h:h
    row = k + 1 + harmonic;
    Xt = Xt + harmonic_matrix(row,:);
end
figure('position', [0, 0, 750, 250]);
plot(time,Xt);
axis tight;
xlabel('Time (sec)');
ylabel('Amplitude');
title(['Synthesized Square Wave for ', num2str(h),' Harmonics Sampled at ', num2str(Fs), ' Hz: #3g']);