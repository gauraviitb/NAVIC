%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT: File containing baseband I and Q samples (4 bytes each) of raw   
% signal, as recorded using GNURadio     
% OPERATION: This program performs acqusition for specified satellites and
% tracking for satellites detected from among them 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 % Sampling Frequency (Fs) used to recorded the data
fs        = 4e6;   

% Path to the recorded data file
fid       = fopen('C:\Everything\Mtech IITB\MTP\SignalProcessing_GPS_IRNSS\Data\GPS\6thMay_GPS_4msps_3_good_9_7_23_27_16.dat');
%fid       = fopen('/home/chirag/Documents/DDP2016/Simulation_codes/gnss-sdr/rec_IITB_19th_sep/signal_source_acquired_gnuradio.dat');

% Satellite Vehicle numbers to be searched for
prn = [10,14,18,32];
N   = size(prn, 2);


% Time (sec) offset from beggining of the file, to start processing
T_off     = 6;  

% Intermediate frequency (Fi), if any. Also useful to compensate for
% LO offset greater than the acquisition search range (of 10 kHz)
fi        = 0;  


code_phase_aq = zeros(1,N); % Stores acquired code-phase
f_doppler     = zeros(1,N); % Stores acquired Doppler Shifts
peak_value    = zeros(1,N); % Stores correlation peak values
sat_detect    = zeros(1,N); % Indicates if the satellite is detected or not

T = 1;                      % Coherent Integration time Tint in ms

%% Acqusition
for q =1:N
[sat_detect(q), code_phase_aq(q), f_doppler(q), acquisition_mat, phase, peak_value(q)] = acquisition(fs, fi, prn(q), T, fid, T_off );
end

%% Tracking

if nnz(sat_detect) > 0       % If any satellite is detected
    
    % Proceed only with satellites successfully acquired
    N             = nnz(sat_detect);
    prn           = prn(find(sat_detect>0));
    peak_value    = peak_value(find(sat_detect>0));
    fc_lo         = f_doppler(find(sat_detect>0));
    code_phase_aq = code_phase_aq(find(sat_detect>0));
    
    % Stores code-phase changes made to acquired code phase by the DLL. 
    code_phase    = zeros(1,N);  
    % Stores LO phase
    phase_lo      = zeros(1,N);
    % Coherent Integration time Tint in ms
    T_int         = 1e-3;      
    
    % Generate PRN codes, and sacel them to {+1,-1} from {1,0}
    ca_code    = zeros(N, fs*1e-3);
    for q =1:N
        ca_code(q,:)    = 2*cacode(prn(q),fs/1.023e6)-1;
    end
    
    % Length of Time (sec) of data processed 
    T_sec      = 10;   % Shouldn't be more than size of the captured data
    time_ms    = 0;    % ms time counter variable
    
    lock_measure = zeros(N,1);
    rec_lock     = zeros(N,T_sec/T_int); % Stores lock measure over time
    rec_bits     = zeros(N,T_sec/T_int); % Stores I component out of PLL
    rec_code_p   = zeros(N,T_sec/T_int); % Stores code-phase out of DLL
    lock         = zeros(1,N);           % Has PLL locked? indicator
     % Stores index of first bit boundary. Initialised to -1, to detect
     % if the bit boundary is not found 
    bit_boundary = -1*ones(1,N);
    
    status      = fseek(fid, fs*T_off*2*4, 'bof');
    % Reads the file and separates the I and Q samples into two columns 
    rec_sig_old = fread(fid, [2,fs*T_int], 'float32')';
    
    while time_ms <= T_sec*1000
        
        rec_sig = fread(fid, [2,fs*T_int], 'float32')'; 
        % 2ms data is read with 1ms overlap with previous data
        sig     = [ rec_sig_old ;rec_sig]; 
        
        for q = 1:N
            rec_I = sig(code_phase_aq(q):code_phase_aq(q)+fs*T_int-1,1);
            rec_Q = sig(code_phase_aq(q):code_phase_aq(q)+fs*T_int-1,2);
            
            % Early-Late codes. Dither = 2 for fs=4MHz
            code_lo   = circshift(repmat(ca_code(q,:),1,T_int/1e-3),ceil(code_phase(q)),2);
            code_lo_e = circshift(repmat(ca_code(q,:),1,T_int/1e-3),ceil(code_phase(q))-fs/2e6,2);
            code_lo_l = circshift(repmat(ca_code(q,:),1,T_int/1e-3),ceil(code_phase(q))+fs/2e6,2);
            
            % Local Oscillator
            carr_lo = cos(2*pi*(0:1/fs:T_int-(1/fs))*fc_lo(q)+ phase_lo(q)) + 1j*sin(2*pi*(0:1/fs:T_int-(1/fs))*fc_lo(q)+ phase_lo(q));
            
            I_complex = dot(carr_lo.*code_lo, rec_I + 1j*rec_Q );
            I         = real(I_complex);
            Q         = imag(I_complex);
            I_e       = dot(carr_lo.*code_lo_e, rec_I + 1j*rec_Q );
            I_l       = dot(carr_lo.*code_lo_l, rec_I + 1j*rec_Q );
            
            %%%%%%%%%%%%%%%%%%%%%%% DLL Code Correction%%%%%%%%%%%%%%%%%%%
            code_correction = abs(I_l)^2 - abs(I_e)^2;
            code_phase(q)   = code_phase(q) + code_correction/(2*peak_value(q));
            %code_phase(q)   = code_phase(q) - 2*(mod(code_phase(q),fs*1e-3) - mod(code_phase(q), fs*1e-3/2)); %To wrap code_phase into [-2000, 2000] instead of [0 4000]
            
            %%%%%%%%%%%%%%%%%%%%%%% Phase correction %%%%%%%%%%%%%%%%%%%%%%
            phase_lo(q)    = mod(phase_lo(q) + 2*pi*T_int*fc_lo(q), 2*pi);
            % atan discriminator capped at pi/4
            pll_correction = ( abs(atan(Q/I))>pi/4 ) * sign(atan(Q/I))*pi/4 + ( abs(atan(Q/I))<pi/4 ) * atan(Q/I);
            phase_lo(q)    = phase_lo(q) + pll_correction ;   
            fc_lo(q)       = fc_lo(q)    + 2*pll_correction/(pi);
            
            lock_measure(q) = lock_measure(q)*0.9 + (I^2 - Q^2)/(I^2 + Q^2);
            
            time_unit = floor(time_ms/(T_int*1000))+1;
            
            rec_bits(q,time_unit)     = I;
            rec_code_p(q,time_unit)   = code_phase(q);
            rec_lock(q,time_unit)     = lock_measure(q);
            
            % Detecting PLL lock
            if lock(q)==0 && lock_measure(q)>5 && abs(sum(sign(rec_bits(q,max(time_unit-19,1):time_unit)))) >=20
                lock(q) = 1;
                JJ=['Satellite ', prn(q), ' Locked!'];disp(JJ);
            end
            
            % Detecting Bit Boundaries
            if lock(q)==1 && bit_boundary(q)<0 && I*rec_bits(q,time_unit-1) < 0 
                bit_boundary(q) = time_unit; 
            end
            
        end
        
        time_ms = time_ms + T_int*1000;
        rec_sig_old  = rec_sig;
        
    end
    fclose(fid);
    JJ=['Lock Status: ', lock];disp(JJ);
    
%% NAV BIT PROCESSING
    % Choose a reference satellite
    [begin_t , sat_ref] = max(bit_boundary);
    
    if begin_t>=0

        bit_boundary_relative = mod(20 - mod(begin_t - bit_boundary, 20),20);
        
        bits_raw = sign(rec_bits);
        
        % Parity Check Matrix
        PB1 = hexToBin('3B1F3480');
        PB2 = hexToBin('1D8F9A40');
        PB3 = hexToBin('2EC7CD00');
        PB4 = hexToBin('1763E680');
        PB5 = hexToBin('2BB1F340');
        PB6 = hexToBin('0B7A89C0');
        
        % Store preamble positions, every 6 sec
        preamble_pos = zeros(N, ceil(T_sec/6)-1); 
        
        for q = 1:N
            bits = bits_raw(q,  begin_t+bit_boundary_relative(q):begin_t + bit_boundary_relative(q)-1+(T_sec-2)*1000);
            
            % Extract NAV bits from PLL output
            temp   = conv(bits(1:(T_sec-2)*1000),ones(1,20));
            bits   = sign(temp(19:20:(T_sec-2)*1000));
            %bits   = sign(bits(1:20:(T_sec-2)*1000));
            
            preamble = [1 -1 -1 -1 1 -1 1 1]; % GPS Preamble
            % Detect Preamble matches
            dtect_pre= find(abs(conv(bits,flip(preamble)))>=8)-7;
            
            %%%%%%%%%%%%%%%%%%%%%%% Check Parity%%%%%%%%%%%%%%%%%%%%%%%%%
            preambles=[];
            pream_count = 0;
            for p=1:size(dtect_pre,2)
                if(dtect_pre(p)>2 && dtect_pre(p)+29<size(bits,2) )
                    d29      = bits(dtect_pre(p)-2)>0;
                    d30      = bits(dtect_pre(p)-1)>0;
                    word_chk = [0 0 (bits(dtect_pre(p):dtect_pre(p)+29)>0)];
                    parityBits = word_chk(27:32); pb = zeros(1,6);
                    pb(1)        = xor(d29,rem(sum(word_chk & PB1),2));
                    pb(2)        = xor(d30,rem(sum(word_chk & PB2),2));
                    pb(3)        = xor(d29,rem(sum(word_chk & PB3),2));
                    pb(4)        = xor(d30,rem(sum(word_chk & PB4),2));
                    pb(5)        = xor(0,rem(sum(word_chk & PB5),2));
                    pb(6)        = xor(xor(d29,d30),rem(sum(word_chk & PB6),2));
                    % To detect false preambles even after parity check,
                    % check if they are spaced 6 sec apart
                    if(pb==parityBits)
                        if(pream_count>=2)
                            d1 = dtect_pre(p) - preambles(pream_count);      
                            d2 = dtect_pre(p) - preambles(pream_count-1);   
                            if mod(d1,300)+mod(d2,300)==0
                                preambles = [preambles dtect_pre(p)];
                                pream_count = pream_count + 1;
                                disp('Hurray');
                            elseif mod(d1,300)==0
                                 preambles(pream_count-1) = preambles(pream_count);
                                 preambles(pream_count) = dtect_pre(p);
                            elseif mod(d2,300)==0
                                preambles(pream_count) = dtect_pre(p);
                            end
                        else
                            preambles = [preambles dtect_pre(p)];
                            pream_count = pream_count + 1;
                            disp('Hurray'); 
                        end
                    end
                end
                
            end
            % (preamable begins at the 'preambles'th position. 
            % So 'preamble-1' 20ms units must be added to reach preamble) 
            preamble_pos(q,1:size(preambles,2)) = preambles-1; 
        end
        
%         %*******************To Store the Subframes into a csv file ******
%         if (size(preambles)>0)
%             packets = bits(preambles(1):preambles(size(preambles,2))-1);
%             packets = (reshape(packets,[300,size(packets,2)/300])')>0;
%         end
%         figure(); plot(rec_bits);
%         dlmwrite('packets_sat.csv',(packets<0.5)');
%         %****************************************************************
        
        preamble_pos_relative = preamble_pos - repmat(preamble_pos(sat_ref,:),N,1);
        
        L = size(preamble_pos,2);
        
        %%%%%%%%%%%%%%%%%%%%%%%% PsuedoRange Computation %%%%%%%%%%%%%%%%%
        psuedo_ranges = zeros(L,N);
        cp            = zeros(L,N);
        
        for q = 1:N
            for l=1:L
                cp(l,q) = round(code_phase_aq(q) + rec_code_p(q,1+begin_t+bit_boundary_relative(q)+20*preamble_pos(q,l)));
                psuedo_ranges(l,q) = ( preamble_pos_relative(q,l)*20e-3*fs + ( bit_boundary_relative(q))*1e-3*fs + cp(l,q) )/fs;
            end 
        end
    else
        disp('Bit boundary not Found')
    end
    
    disp('Satellites Found are:'); prn
    
    sat_PR = ( 65e-3  +psuedo_ranges)*299792458;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

% 
%             for p=1:size(dtect_pre,2)
%                 if(dtect_pre(p)>2 && dtect_pre(p)+29<size(bits,2) )
%                     d29      = bits(dtect_pre(p)-2);
%                     d30      = bits(dtect_pre(p)-1);
%                     word_chk = [d29 d30 (bits(dtect_pre(p):dtect_pre(p)+29))];
%                     pb = mod(H*word_chk',4)-2;
%                     parityBits = word_chk(27:32);
%                     if(pb'==parityBits)
%                         if(pream_count>=2)
%                             d1 = dtect_pre(p) - preambles(pream_count);      % To take care of false preamble detection even after FEC
%                             d2 = dtect_pre(p) - preambles(pream_count-1);   
%                             if mod(d1,300)+mod(d2,300)==0
%                                 preambles = [preambles dtect_pre(p)];
%                                 pream_count = pream_count + 1;
%                                 disp('Hurray');
%                             elseif mod(d1,300)==0
%                                  preambles(pream_count-1) = preambles(pream_count);
%                                  preambles(pream_count) = dtect_pre(p);
%                             elseif mod(d2,300)==0
%                                 preambles(pream_count) = dtect_pre(p);
%                             end
%                         else
%                             preambles = [preambles dtect_pre(p)];
%                             pream_count = pream_count + 1;
%                             disp('Hurray'); 
%                         end
%                     end
%                 end
%                 
%             end