fs        = 4e6;    % Sampling frequency 
% ?? What's the right sampling frequency to choose? As fs is increased the SNR
% reduces, but does it aid pseudo code de-corelation ??

T         = 2 ;     % length of signal analysed for acqusition in ms (should be integer) % ** T should be < 1/(4*freq_bin)**
% For a given fs, it seems increasing T brings advantage. However, the freq-bin needs to be reduced
% accordingly. Thus it's good to begin with T = 1 and then move to higher T
% in Tracking. ?? Increasing T for better detection probablity vs more
% expensive acquition (in terms of time, power). How to settle tradeoff ??
                   
%fid       = fopen('/home/chirag/Documents/DDP2016/Simulation_codes/Matlab_C++_Progress/September/Acquisition/data/gps_data_recIITB_3rd_sep.dat');
%fid = fopen('/media/chirag/1CC0-29C2/3.dat');
fid       = fopen('/home/chirag/Documents/DDP2016/Simulation_codes/Matlab_C++_Progress/March/data/March7_SAMEER_315PM/irnss_rhcp1_4msps_1.dat');
% [sats 10,8,18,11,32 on 19th Sep 19th September, 45.5sec use sat
% 10,11,32,14], See T_off = 45 and prn = 18 to see effects of having nearly zero doppler
% 3rd sep, T_off = 5, prn=10, repeated acqusition helps in locking the
% satellite inspite of bad acqusition peak

T_off     = 7; % Time offset in seconds

fi  = 0;  % Intermediate freq at which signal is received from frontend
% ?? Advantages of processing at IF vs at Baseband ?? Baseband processing
% requires I and Q values. What's the advantage of using both I and Q
% values of the sampled signal ??

prn = [3];
no_sats = size(prn,2);

sat_detect = zeros(1,no_sats);
code_phase_aq = zeros(1,no_sats);
 f_doppler = zeros(1,no_sats);

for q=1:no_sats
[sat_detect(q), code_phase_aq(q), f_doppler(q) ] = acquisition(fs, fi, prn(q), T, fid, T_off );
end
% ?? How to detect lock in acquisition??
%%
if sat_detect >0

code_phase = code_phase_aq;
fc_lo      = f_doppler; 
phase_lo   = 0;
T_int      = 1e-3;  %Integration time to find ip and qp
T          = 1000; % length of signal processed in ms
ca_code    = 2*cacode_irnss(prn,fs/1.023e6,'s')-1;
rec_bits   = [];
rec_freq   = [];
rec_code_p = [];
rec_peaks  = [];
rec_I_comp = [];
rec_freq_er= [];
rec_code_er= [];
preambles  = [];
lock_measure = 10;%pi
lock_or_not  = 0;

T_sec      = 20;   % Shouldn't be more than the size of the captured data
for sec = 1:T_sec
    sec
    %if sec>5                    % Adaptive increase of corelation length
    %    T_int = 2e-3;
    %end
    rec_sig = fread(fid, [2,fs*1e-3*T], 'float32')';
    bpsk_I = rec_sig(:,1);
    bpsk_Q = rec_sig(:,2);


    code_lo   = circshift(repmat(ca_code,1,T_int/1e-3),[0,ceil(code_phase)]);
    code_lo_e = circshift(repmat(ca_code,1,T_int/1e-3),[0,ceil(code_phase)-2]);
    code_lo_l = circshift(repmat(ca_code,1,T_int/1e-3),[0,ceil(code_phase)+2]);
    %**** PLL*****

    code_rate = 1023*T_int*1e3 + ((fs/1.023e6)*fc_lo/1575.42)*T_int ; %per ms

    N     = floor(T/20);

    I            = zeros(N*20*1e-3/T_int,1);
    Q            = zeros(N*20*1e-3/T_int,1);
    phase_record = zeros(N*20*1e-3/T_int,1);
    fc_lo_record = zeros(N*20*1e-3/T_int,1);
    code_phi_rec = zeros(N*20*1e-3/T_int,1);
    peak_record  = zeros(N*20*1e-3/T_int,1);
    error        = zeros(N*20*1e-3/T_int,1);
    code_err     = zeros(N*20*1e-3/T_int,1);

    for p = 1:N*20/(T_int*1e3)   % every ms
        
        carr_lo = cos(2*pi*(0:1/fs:T_int-(1/fs))*fc_lo+ phase_lo).*code_lo + 1j*sin(2*pi*(0:1/fs:T_int-(1/fs))*fc_lo+ phase_lo).*code_lo;
        
        carr_lo_e = cos(2*pi*(0:1/fs:T_int-(1/fs))*fc_lo+ phase_lo).*code_lo_e + 1j*sin(2*pi*(0:1/fs:T_int-(1/fs))*fc_lo+ phase_lo).*code_lo_e;
        
        carr_lo_l = cos(2*pi*(0:1/fs:T_int-(1/fs))*fc_lo+ phase_lo).*code_lo_l + 1j*sin(2*pi*(0:1/fs:T_int-(1/fs))*fc_lo+ phase_lo).*code_lo_l;
    
        I_complex = dot(carr_lo, bpsk_I( (p-1)*fs*(T_int)+1 : p*fs*(T_int)) + 1j*bpsk_Q( (p-1)*fs*(T_int)+1 : p*fs*(T_int)) ); % integration to get ip, qp
        I(p) = real(I_complex); Q(p)= imag(I_complex);
        I_e  = dot(carr_lo_e, bpsk_I( (p-1)*fs*(T_int)+1 : p*fs*(T_int)) + 1j*bpsk_Q( (p-1)*fs*(T_int)+1 : p*fs*(T_int)) );
        I_l  = dot(carr_lo_l, bpsk_I( (p-1)*fs*(T_int)+1 : p*fs*(T_int)) + 1j*bpsk_Q( (p-1)*fs*(T_int)+1 : p*fs*(T_int)) );
        
        code_correction = abs(I_l)^2 - abs(I_e)^2;
        code_err(p)      =  code_correction;
        phase_lo = mod(phase_lo + 2*pi*T_int*fc_lo, 2*pi); 
        
        % PLL correction
        % Cap the correction so that large noise doesn't flip the tracker
        % to theta+pi (since atan is agnostic to pi difference). Such shift
        % cause random flipping of received nav bits
        correction = ( abs(atan(I(p)/Q(p)))>pi/4 ) * sign(atan(I(p)/Q(p)))*pi/4 + ( abs(atan(I(p)/Q(p)))<pi/4 ) * atan(I(p)/Q(p));
        %correction = atan(I(p)/Q(p));
        
        %lock_measure = lock_measure*0.95 + correction;
        lock_measure = lock_measure*0.9 + abs(I(p)/Q(p));
        if lock_or_not == 0 && lock_measure <10
            disp('Locked!')
            lock_or_not = 1;
        end
        
        if lock_or_not == 1 && lock_measure >10
            disp('NOT Locked!')
            lock_or_not = 0;
        end
        %error(p) = correction;
        error(p) = lock_measure; 
        
        phase_lo = phase_lo - correction ;    % atan discriminator
        fc_lo    = fc_lo - 2*correction/(pi);

        
        code_phase = mod(code_phase + code_correction/(8000*(T_int/1e-3)^2), fs*1e-3);   % 4000 is choosen as half of max peak value
        code_lo   = circshift(repmat(ca_code,1,T_int/1e-3),[0,floor(code_phase)]);
        code_lo_e = circshift(repmat(ca_code,1,T_int/1e-3),[0,floor(code_phase)-2]);
        code_lo_l = circshift(repmat(ca_code,1,T_int/1e-3),[0,floor(code_phase)+2]);

        code_phi_rec(p) = floor(code_phase);
        phase_record(p) = phase_lo;
        fc_lo_record(p) = fc_lo;
        
    end
    
rec_bits = [rec_bits Q'];
rec_I_comp = [rec_I_comp I'];
rec_code_p = [rec_code_p code_phi_rec'];
rec_freq   = [rec_freq fc_lo_record'];
p          = mean(I.^2 + Q.^2);
rec_peaks = [rec_peaks p ];
rec_freq_er = [rec_freq_er, error'];
rec_code_er = [rec_code_er code_err'];
% pause
% plot(rec_bits)

end
fclose(fid);
%% Nav Bit processing

bits_raw = sign(rec_bits);
%plot(bits_raw)

figure(); plot(rec_bits);
 figure();plot(rec_freq_er);grid on
end
