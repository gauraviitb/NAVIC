clc
clear
fs = 1;
fi = 0;
T = 1;
a = cacode(7, fs);    % PRN sequence from SV 7
incoming = cacode(7, fs); 
incoming = circshift(incoming', 421) ;
incoming = -1.*cos(2*pi*(fi+5000)*T).*incoming';  % incoming data from SV 7
%a = [1 2 3 4 5 6];  % PRN
%incoming = [6 5 4 3 2 1]; % incoming data

freq_bin   = 500/T; 
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
max_index_fd = index1(index2) 
%max_fd = fd_lo(max_index_fd) ;
%max_t = max_index_t/fs ;
max_ele = acq_mat_ele(max_index_fd, max_index_t);

f_doppler              = -freq_range + freq_bin*(max_index_fd-1)
code_phase_aq          = (0 + code_bin*(max_index_t-1))

if(max_ele < 100000 && max_ele > -100000)
disp('Satellite not found. Signal not good.');
else
        if (max_ele < 0)
          bit = 0
        else bit = 1
        end;
 disp('Acqusition Complete. Satellite found');
 end