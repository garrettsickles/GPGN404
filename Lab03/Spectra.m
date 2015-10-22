function values = Spectra(frequency, value, type)
%Spectra - Plots phase and amplitude

values = zeros(2, length(frequency)*2-1);
domain = zeros(1, length(frequency)*2-1);
range = zeros(1, length(frequency)*2-1);

j = 1;
for i = 1:length(frequency)
    if frequency(i) >= 0
        domain(j) = frequency(i);
        if strcmp(type, 'Amplitude')
            if frequency(i) > 0
                domain(j+1) = (-1.0)*domain(j);
                range(j) = abs(value(i));
                range(j+1) = abs(value(i));
                j = j+1;
            elseif frequency(i) == 0
                range(j) = abs(value(i));
            end
            j=j+1;
        elseif strcmp(type, 'Phase')
            range(j) = atan(imag(value(i))/real(value(i)));
            if frequency(i) > 0
                domain(j+1) = (-1.0)*domain(j);
                range(j) = atan(imag(value(i))/real(value(i)));
                range(j+1) = -atan(imag(value(i))/real(value(i)));
                j = j+1;
            elseif frequency(i) == 0
                range(j) = atan(imag(value(i))/real(value(i)));
            end
            j=j+1;
        end 
    end
end
figure;
stem(domain, range);
xlabel({'Frequency (Hz)'});
title(strcat(type, ' Spectrum'));
ylabel({type});
end

