clc
clear
fs = 1;
fi = 0;

code_phase_shift = 9;
doppler_shift = 800;
data_bit = -1;    %1 or -1
snr = -20;     % SNR per sample in dB


T = linspace(0, 0.001,1023*20);
a = cacode(8, fs);    % PRN sequence at receiver from SV 8
a = cat(2,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a) ;
%a = cat(2,a,a,a,a,a);
incoming = cacode(8, fs); % incoming PRN sequence
incoming = cat(2, incoming,incoming, incoming, incoming, incoming, incoming, incoming, incoming, incoming, incoming, incoming, incoming, incoming, incoming, incoming, incoming, incoming, incoming, incoming, incoming);
%incoming = cat(2, incoming,incoming, incoming, incoming, incoming);
incoming = circshift(incoming', code_phase_shift) ;
     
carrier = cos(2*pi*(fi + doppler_shift).*T);    % doppler freq 
subplot(221);
plot(T,incoming);
subplot(222);
plot(T,carrier);
hold on;
carrier = sign(carrier);
plot(T,carrier);
hold off;
%incoming = circshift(incoming', 50) ;       % incoming delay
incoming = data_bit*carrier.*incoming';
subplot(223);
plot(T, incoming);



y = awgn(incoming,snr);  % incoming data from SV 7 with noise
incoming = y;
subplot(224);
plot(T,incoming);


freq_bin   = 500;
freq_range = 10e3;
code_bin    = fs/1e6;
fd_lo   = -freq_range:freq_bin:freq_range;
b = a;
for j = 1 : size(a,2)
b = shift(b);  % circular shift by 1
for i = 1 : size(fd_lo,2)
[d , c] = scale (fi , fd_lo(i) , T , b );
acq_mat_ele(i,j) = correlate (d, c, incoming);
end
i = 1;
end

[max_ele, index1] = max(abs(acq_mat_ele));
[max_ele, index2] = max(abs(max_ele)) ;
max_index_t = index2
max_index_fd = index1(index2);
max_ele = acq_mat_ele(max_index_fd, max_index_t);

f_doppler              = -freq_range + freq_bin*(max_index_fd-1)
code_phase_aq          = (0 + code_bin*(max_index_t-1))

if(max_ele < 10000000 && max_ele > -10000000)
disp('Satellite not found. Signal not good.');
else
        if (max_ele < 0)
          bit = 0
        else bit = 1
        end;
 disp('Acqusition Complete. Satellite found');
 end
