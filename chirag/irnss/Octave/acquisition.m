function [sat_detect, code_phase_aq, f_doppler] = acquisition(fs, fi, prn, T, fid, T_off )

if T_off > 0
    status = fseek(fid, fs*T_off*2*4, 'bof'); % Offset by T sec from beginning of file (each sample is 4 bytes real and 4 bytes imaginary)
    if status ~= 0
     disp('Change T_off so that 8*fs*T_off is an (even) integer');
     return
    end
end

freq_bin   = 500/T; % As T is increased freq_bin should reduce 
freq_range = 10e3; 

fd_lo   = -freq_range:freq_bin:freq_range; 
M       = size(fd_lo,2);    % M = no of freq bin to search 

code_lo = 2*repmat(cacode_irnss(prn,fs/1.023e6,'s'),1,T)-1;
N       = size(code_lo,2);  % N = sequence length per freq

code_bin = fs/1e6;           
code_phi_lo = 0:code_bin:N/T;
L           = size(code_phi_lo,2); % L = Number of code phases to search

peak_threshold  = 10000; % should be changed as fs changes or input signal power changes
% ?? Decide peak threshold for a given false alarm proablity (pg 254 GPS: Theroy and APplications, volume 1)

peak_value      = 0;
ms_count        = 0; % ms count until satellite is detected
lock            = 0;
max_index_old       = 0;
%% Acquisition
    rec_sig   = fread(fid, [2,fs*1e-3*T], 'float32')';
    rec_sig_I = rec_sig(:,1)';
    rec_sig_Q = rec_sig(:,2)';    

    I_lo = cos(2*pi*(fi + fd_lo)'*(1/fs:1/fs:T*1e-3));
    Q_lo = sin(2*pi*(fi + fd_lo)'*(1/fs:1/fs:T*1e-3));

while lock==0 && size(rec_sig,1)== fs*T*1e-3 &&  ms_count<10

    ms_count = ms_count + 1
    
    acquisition_mat = zeros(M, L); % Search space

    for v = 1:L
        temp_I = I_lo.*repmat( circshift(code_lo,[0,code_phi_lo(v)]),M,1);
        temp_Q = Q_lo.*repmat( circshift(code_lo,[0,code_phi_lo(v)]),M,1);
        I      = temp_I*rec_sig_I' + temp_Q*rec_sig_Q';
        Q      = -temp_Q*rec_sig_I' + temp_I*rec_sig_Q';
        acquisition_mat(:,v)= I.^2 + Q.^2;
    end

    [peak_value, max_index] = max(acquisition_mat(:))
    [f_d,code_index]       = ind2sub(size(acquisition_mat),max_index);
    rec_sig   = fread(fid, [2,fs*1e-3*T], 'float32')';
    rec_sig_I = rec_sig(:,1)';
    rec_sig_Q = rec_sig(:,2)';
    
    if max_index == max_index_old || peak_value >= peak_threshold 
        lock = 1;
    end
    
    max_index_old = max_index;
    
    figure()
    plot(1:M*L, reshape(acquisition_mat/T,[1,M*L]))
end

if  lock==1
    sat_detect = ms_count*T;
    f_doppler              = -freq_range + freq_bin*(f_d-1);
    code_phase_aq          = (0 + code_bin*(code_index-1));
    disp('Acqusition Complete. Satellite found')
else
    sat_detect = 0;
    f_doppler  = 0;
    code_phase_aq = 0;
    disp('Satellite not found. Signal not good.')
end 

%fclose(fid);
%% Plot
% figure()
% surf(acquisition_mat/T,gradient(acquisition_mat))
% rotate3d on
peak_value
figure()
plot(1:M*L, reshape(acquisition_mat/T,[1,M*L]))
acquisition_mat=reshape(acquisition_mat/T,[1,M*L]);
end
