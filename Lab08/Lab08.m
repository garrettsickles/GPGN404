function Lab08()
    Problem1();
end

function Problem1()
    Fs = 25.0;
    Ts = 1/Fs;
    
    Tmin = 0.0;
    Tmax = 0.5;
    
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
    
    figure;
    plot(tc, y, 'color', 'b');
    hold on;
    plot(tc, z, 'color', [1,0.3,0]);
    hold on;
    plot(ts, x, 'o', 'color', 'k');
    xlabel('Time (sec)');
    ylabel('Amplitude');
    legend('Original','Reconstructed','Samples Points');
    title('Problem 1');
end

% function Problem2()
%     Fs = 15.0;
%     Ts = 1/Fs;
%     
%     Tmin = 0.0;
%     Tmax = 0.5;
%     t = Tmin:Ts:Tmax;
%     n = Tmin/Ts:1:Tmax/Ts;
%     N = length(n);
%     tc = Tmin:0.001:Tmax;
%     
%     x = (2.0)*cos(3.1415926*(10*n*Ts+1/6))+cos(3.1415926*(20*n*Ts+1/3));  
%     y = (2.0)*cos(3.1415926*(10*tc+1/6))+cos(3.1415926*(20*tc+1/3));
%     z = zeros(length(tc),1);
%     
%     for i=1:length(tc)
%         for j=1:N
%             z(i) = z(i) + x(j) * sinc((tc(i)-t(j))/Ts);
%         end
%     end
%     
%     figure;
%     plot(tc, y, 'color', 'b');
%     hold on;
%     plot(tc, z, 'color', [1,0.3,0]);
%     hold on;
%     plot(Ts.*n, x, 'o', 'color', 'k');
%     xlabel('Time (sec)');
%     ylabel('Amplitude');
%     legend('Original','Reconstructed','Samples Points');
%     title('Problem 2');
% end