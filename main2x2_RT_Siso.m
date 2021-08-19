
%%   MAIN FUNCTION
function main2x2_RT(app2,file)

%% Load Data
app2.enable.HighpassEqui = 1;   % remove high freq comp from raw data
[complex_data,app2] = load_data(app2, file);

%% controls / select

    app2.Sweeps2Fetch = 16;
    app2.Sweeps_inFrame = 64;
    app2.Nfftr = app2.nSamp;
    app2.Nfftv = 2*app2.Sweeps_inFrame;

    app2.HanningWindowDopp = 0;
    app2.HanningWindowRange = 1;
    app2.RemoveClutter = 0;
    
    app2.RangeNotchFilter = 1; % range filter
    app2.RangeNotch = 4000;

    app2.DopplerNotchFilterCheckBox.Value = 0; % Doppler filter    

%% range profile segments

    SeeAllRanges = 0;   % 0 = seperate ranges; 1 = combined ranges
    
    % Tx2-7m, Rx2-7m cable
    if(app2.bw==1e9)
        app2.R1_start_m = 0;
        app2.R1_stop_m = 50;
        %app2.R2_start_m = 11.5;%3;
        %app2.R2_stop_m = 9.6*2;%4;
        %interference link
        %app2.R3_start_m = 5.85;
        %app2.R3_stop_m = 11.5;
    elseif(app2.bw==250e6)
        app2.R1_start_m = 0.6;0.6;
        app2.R1_stop_m = 50;40;
        %app2.R2_start_m = 12.6;%3;
        %app2.R2_stop_m = 9.6*2;%4;
        %interference link
        %app2.R3_start_m = 7.2;
        %app2.R3_stop_m = 12.6;    
    end

    if(SeeAllRanges)
        % all ranges, to show interference, or for analysis
        app2.R1_stop_m = 50; app2.R1_start_m = 3e8/2/app2.bw;
        %app2.R2_stop_m = 9.6*2; app2.R2_start_m = 3e8/2/app2.bw;
        %app2.R3_stop_m = 9.6*2; app2.R3_start_m = 3e8/2/app2.bw;
        
        % Range Adjustment for x & y radars
        app2.R1_rangeOffset = 0; % 1.88->xmax60,2x2
        %app2.R2_rangeOffset = 0;
        %app2.R3_rangeOffset = 0;  
    else
        % Range Adjustment for x & y radars
        app2.chAstartRange = 0; 1;  % fixed offset
        app2.chBstartRange = 0; 7;  % fixed offset
        app2.chCstartRange = 0; 4;  % fixed offset
        app2.R1_rangeOffset = app2.R1_start_m - app2.chAstartRange;1.188; % 1.88->xmax60,2x2
%        app2.R2_rangeOffset = app2.R2_start_m - app2.chBstartRange;
%        app2.R3_rangeOffset = app2.R3_start_m - app2.chCstartRange;
    end
    
    % Doppler Axis
    app2.toDoppAxis = 2*app2.fc/3e8;

%% graphs
    
    min_d = 35;
    max_d = 100;
    commonGaindB = 0;
    %graph 1
    app2.time_offs = 0;
    app2.dBGain1 = 1 + commonGaindB;   % dB gain adjustment
    app2.clim_min_rd1 = 15;
    app2.clim_max_rd1 = 55;%% range doppler 
    app2.clim_min_d1 = min_d;
    app2.clim_max_d1 = max_d;%% doppler   
    app2.clim_min_r1 = 35;41;44;10;
    app2.clim_max_r1 = 61;70;%% range
    %graph 2
    app2.dBGain2 = 12 + commonGaindB;   % dB gain adjustment
    app2.clim_min_rd2 = 0;
    app2.clim_max_rd2 = 55;%% range doppler 
    app2.clim_min_d2 = 28;min_d-2;
    app2.clim_max_d2 = 48;max_d-12;%% doppler   
    app2.clim_min_r2 = 25;30;43;10;
    app2.clim_max_r2 = 50;58;%% range
    %graph 3 interferance h12,h21
    app2.dBGain3 = 0 + commonGaindB;   % dB gain adjustment
    app2.clim_min_rd3 = 15;
    app2.clim_max_rd3 = 55;%% range doppler 
    app2.clim_min_d3 = min_d;
    app2.clim_max_d3 = max_d;%% doppler   
    app2.clim_min_r3 = 35;30;43;10;
    app2.clim_max_r3 = 61;58;%% range

    app2.doppScFactMeas = 1; %Doppler-graph settings for paper
    app2.time_offset = 18;47.5;      % 18: xmax ; 47.5: elliptical
    app2.dopp_scalingAnal = 1;0.84;

    app2.gr_color = 50*[1 1 1]/255;
    app2.DSindx = 7;     % downSampleIndex
    app2.gr_width = 0.23;0.715/2;    % OfficeMonitor = 0.53(for subplot); 0.23(singlePlot) , Laptop = 0.715/2;
    app2.gr_len = 0.3;0.41;       % OfficeMonitor = 0.3 , Laptop = 0.41;
    app2.ltx_fontSz = 12;
    
    % graph enables
    app2.enable.R=1;    %range
    app2.enable.D=1;    %doppler

%% while loop
complex_data_Stitch = [];
ii_data = 1;

while(1)
    time_v = [1:ii_data] * app2.sTime * app2.Sweeps2Fetch .* app2.nTx;
    time_v= time_v - app2.time_offs;
    
    %% data, sliding window
    if(app2.Sweeps2Fetch <= app2.Sweeps_inFrame)
        if(ii_data==1)  % Initial State / initial data 256 sweeps
            complex_data_a = complex_data( (ii_data-1) * app2.Sweeps2Fetch * app2.nSamp + 1 : ii_data * app2.Sweeps2Fetch * app2.nSamp , :);
            ii_data = 2;
            angle_v_t=[];
            Rng_psd_dB = [];
            Dopp_psd_dB = [];
            Rng_psd_dB_2 = [];
            Dopp_psd_dB_2 = [];
            Rng_psd_dB_3 = [];
            Dopp_psd_dB3_ = [];
            Rng_psd_dB_4 = [];
            Dopp_psd_dB_4 = [];
        elseif (ii_data * app2.Sweeps2Fetch * app2.nSamp < length(complex_data)) % send 256 sweeps down in each iteration
            complex_data_a = complex_data( (ii_data-1) * app2.Sweeps2Fetch * app2.nSamp + 1 : ii_data * app2.Sweeps2Fetch * app2.nSamp , :);
            ii_data = ii_data + 1;
        else
            app2.PauseAppCheckBox.Value = 1;
            disp('end of file reached');
            break;
            angle_v_t=[];
            Rng_psd_dB = [];
            Dopp_psd_dB = [];
            Rng_psd_dB_2 = [];
            Dopp_psd_dB_2 = [];
            Rng_psd_dB_3 = [];
            Dopp_psd_dB3_ = [];
            Rng_psd_dB_4 = [];
            Dopp_psd_dB_4 = [];
        end    
    end
    
    if(size(complex_data_Stitch,1) < app2.Sweeps_inFrame * app2.nSamp)
        complex_data_Stitch = [complex_data_Stitch; complex_data_a];
    else
        complex_data_Stitch = [complex_data_Stitch(app2.Sweeps2Fetch * app2.nSamp + 1 : end,:) ; complex_data_a];
    end                                
    df_rng = app2.sampRate /  app2.Nfftr; % frequency resolution in range
    df_dopp = app2.prf / app2.Nfftv; % frequency resolution in Doppler

    nChirps = floor(size(complex_data_Stitch,1)/app2.nSamp);  
    
    %% for loop 
    for ch_i = 1 : size(complex_data_Stitch,2)      % ch_i is Channel # ; or for dbf->ScanAngle#
        
        complex_data_1 = reshape(complex_data_Stitch(1:app2.nSamp*nChirps,ch_i),app2.nSamp,nChirps);
        
        %% hanning window & clutter
        Nframe = app2.Sweeps_inFrame;
        if app2.HanningWindowDopp
            my_Win_v = hann(size(complex_data_1,2));
            my_Win = transpose(repmat(my_Win_v,[256,1]));
        elseif app2.HanningWindowRange
             my_Win_v = hann(size(complex_data_1,1));
             my_Win = repmat(my_Win_v,[1,Nframe]);
            
            %my_Win_v = chebwin(size(complex_data_1,1),800);
            %my_Win = repmat(my_Win_v,[1,Nframe]);
        end
        
        dopp_axis = linspace(-app2.prf/2,app2.prf/2,app2.Nfftv)*3e8/app2.fc/2;
        rng_axis = linspace(0, app2.sampRate/2, app2.Nfftr/2+1)*3e8* ...
            app2.sTime/(2*app2.bw);

        if app2.HanningWindowRange
            if(size(complex_data_1,2)==size(my_Win,2))
                data_matrix_1 = complex_data_1.*(my_Win);
            else
                data_matrix_1 = complex_data_1;
            end
        else
            data_matrix_1 = complex_data_1;
        end
        if(app2.RemoveClutter)
            data_matrix_1 = bsxfun(@minus, data_matrix_1, mean(data_matrix_1, 2)); % remove static objects  
        end
        
        %% Filtering in range
        if app2.RangeNotchFilter
            notchfreq = app2.RangeNotch+1;
            if notchfreq >= app2.sampRate/2.5
                notchfreq = app2.sampRate/2.5;
            end
            F1 = 2*notchfreq/app2.sampRate;
            F2 = 2*app2.sampRate/2.1/app2.sampRate;
            F12 = [F1,F2];
            N = 8;
            [BF,AF] = butter(N,F12);
            data_matrix_1 = filter(BF,AF,data_matrix_1,[],1);
        end

        %% Range profiles
        Rng_dft_1a = fft(data_matrix_1, app2.Nfftr, 1);
        rng_axis_delta = rng_axis(2);
        if(ch_i == 1)
            app2.start1EditField.Value = round(app2.R1_start_m / rng_axis_delta);
            app2.stop1EditField.Value = round(app2.R1_stop_m / rng_axis_delta);
            Rng_dft_1 = Rng_dft_1a(app2.start1EditField.Value:app2.stop1EditField.Value, :);
            rng_axis1 = rng_axis(1 : app2.stop1EditField.Value-app2.start1EditField.Value + 1);
            
            % h12,h21 link
%            app2.start3EditField.Value = round(app2.R3_start_m / rng_axis_delta);
%            app2.stop3EditField.Value = round(app2.R3_stop_m / rng_axis_delta);
%            Rng_dft_interf = Rng_dft_1a(app2.start3EditField.Value:app2.stop3EditField.Value, :);
%            rng_axis3 = rng_axis(1 : app2.stop3EditField.Value-app2.start3EditField.Value + 1);
            
            rng_axis = rng_axis1;
        elseif(ch_i == 2)
            app2.start2EditField.Value = round(app2.R2_start_m / rng_axis_delta);
            app2.stop2EditField.Value = round(app2.R2_stop_m / rng_axis_delta);
            Rng_dft_1 = Rng_dft_1a(app2.start2EditField.Value:app2.stop2EditField.Value, :);
            rng_axis2 = rng_axis(1 : app2.stop2EditField.Value-app2.start2EditField.Value + 1);
            rng_axis = rng_axis2;
        end

        %% Filtering in Doppler
        if app2.DopplerNotchFilterCheckBox.Value
            notch_vr = app2.DopplerNotchSlider.Value+0.01;
            notch_doppler = 2*notch_vr/(3e8/app2.fc);
            stopfreq = notch_doppler;
            if stopfreq > app2.prf/3
                stopfreq = app2.prf/3.1;
            end
            passfreq = 1.5*stopfreq;
            ws = stopfreq/(app2.prf/2);
            wp = passfreq/(app2.prf/2);
            Rp = 2;
            Rs = 30;
            [N,Wn] = buttord(wp,ws,Rp,Rs);
            [num,den] = butter(N,Wn,'high'); % Highpass filter
            Rng_dft_1 = filter(num,den,Rng_dft_1,[],2); % Doppler Notch filter
        end

        %% range Doppler PSD
        RngDopp_psd_1 = abs(fftshift(fft(Rng_dft_1,app2.Nfftv,2),2)).^2/(df_dopp*df_rng);
%         RngDopp_psd_dB = 10*log10(RngDopp_psd_1);
        % for h12,h21 (interference link)
%        RngDopp_psd_interf = abs(fftshift(fft(Rng_dft_interf,app2.Nfftv,2),2)).^2/(df_dopp*df_rng);
%         RngDopp_psd_dB_interf = 10*log10(RngDopp_psd_interf);
         
        %% Doppler PSD
        if ch_i == 1
            Dopp_psd_dB1(ii_data-1,:) = (sum(RngDopp_psd_1))./sum(sum(RngDopp_psd_1)); %+ app2.dBGain1;
%            Dopp_psd_dB3(ii_data-1,:) = 10*log10(sum(RngDopp_psd_interf)) + app2.dBGain3;
            
        elseif ch_i == 2
            Dopp_psd_dB2(ii_data-1,:) = (sum(RngDopp_psd_1))./sum(sum(RngDopp_psd_1));% + app2.dBGain2;
        end
        
        %% Range PSD
        if ch_i == 1
            Rng_psd_dB1(ii_data-1,:) = 10*log10(sum(transpose(RngDopp_psd_1)));
%            Rng_psd_dB3(ii_data-1,:) = 10*log10(sum(transpose(RngDopp_psd_interf)));
            
        elseif ch_i == 2
            Rng_psd_dB2(ii_data-1,:) = 10*log10(sum(transpose(RngDopp_psd_1)));
        end
        
    end % end of for
%     elp_time(ii_data)=toc;
end % end of while

% avg_whileLoop_time = mean(elp_time)

    %% Doppler PSD
    if(app2.enable.D==1)
        figure(211), set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.02, 0.68-0*app2.gr_len, app2.gr_width, app2.gr_len]);
        app2.clim_min_d1 = 0;
        app2.clim_max_d1 = max(max(Dopp_psd_dB1));
        psd_plot(app2,'doppPSD',time_v(1:ii_data-1), -dopp_axis*app2.doppScFactMeas, transpose(Dopp_psd_dB1),app2.clim_min_d1, app2.clim_max_d1,'Time, $t$ (s)','Radial velocity, $\dot{d}_{ij}(t)$ (m/s)',app2.ltx_fontSz)

        % for interf link
%        figure(212), set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.02, 0.68-1*app2.gr_len, app2.gr_width, app2.gr_len]);
%        psd_plot(app2,'doppPSD',time_v(1:ii_data-1), -dopp_axis*app2.doppScFactMeas, transpose(Dopp_psd_dB3),app2.clim_min_d3, app2.clim_max_d3,'Time, $t$ (s)','Radial velocity, $\dot{d}_{ij}(t)$ (m/s)',app2.ltx_fontSz)

%        figure(222), set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.02, 0.68-2*app2.gr_len, app2.gr_width, app2.gr_len]);
%        psd_plot(app2,'doppPSD',time_v(1:ii_data-1), -dopp_axis*app2.doppScFactMeas, transpose(Dopp_psd_dB2),app2.clim_min_d2, app2.clim_max_d2,'Time, $t$ (s)','Radial velocity, $\dot{d}_{ij}(t)$ (m/s)',app2.ltx_fontSz)    
    end
    
    %% Range PSD
    if(app2.enable.R==1)
        % range axis offset calibration
        rng_axis1 = rng_axis1 + app2.R1_rangeOffset;
%        rng_axis2 = rng_axis2 + app2.R2_rangeOffset;
%        rng_axis3 = rng_axis3 + app2.R3_rangeOffset;

        figure(311), set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5, 0.68-0*app2.gr_len, app2.gr_width, app2.gr_len]);
    %     psd_plot(app,'rngPSD',time_v(1:ii_data-1), rng_axis1(1:end-1), transpose(Rng_psd_dB1),app2.clim_min_r1, app2.clim_max_r1,'Time, $t$ (s)','Range, $d_{ij}(t)$ (m)',app2.ltx_fontSz)
        psd_plot(app2,'rngPSD',time_v(1:ii_data-1), rng_axis1(1:end), transpose(Rng_psd_dB1),app2.clim_min_r1, app2.clim_max_r1,'Time, $t$ (s)','Range, $d_{ij}(t)$ (m)',app2.ltx_fontSz)

        %for interf link h12
%        figure(312), set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5, 0.68-1*app2.gr_len, app2.gr_width, app2.gr_len]);
%        psd_plot(app2,'rngPSD',time_v(1:ii_data-1), rng_axis3(1:end), transpose(Rng_psd_dB3),app2.clim_min_r3, app2.clim_max_r3,'Time, $t$ (s)','Range, $d_{ij}(t)$ (m)',app2.ltx_fontSz)

%        figure(322), set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5, 0.68-2*app2.gr_len, app2.gr_width, app2.gr_len]);
    %     psd_plot(app,'rngPSD',time_v(1:ii_data-1), rng_axis2(1:end-1), transpose(Rng_psd_dB2),app2.clim_min_r2, app2.clim_max_r2,'Time, $t$ (s)','Range, $d_{ij}(t)$ (m)',app2.ltx_fontSz)
%        psd_plot(app2,'rngPSD',time_v(1:ii_data-1), rng_axis2(1:end), transpose(Rng_psd_dB2),app2.clim_min_r2, app2.clim_max_r2,'Time, $t$ (s)','Range, $d_{ij}(t)$ (m)',app2.ltx_fontSz)
    end

end
%% Other Functions
function [complex_data,app2] = load_data(app2, file)
%% load data

radarDat = load(file);

app2.nRx = radarDat.nRx;
app2.nTx = radarDat.nTx;
app2.fStart = radarDat.fStart;
app2.fStop = radarDat.fStop;
app2.sTime = radarDat.sTime;%*1e-3;
app2.nSamp = radarDat.nSamp;
rawdata = radarDat.ADCsamples;
app2.nSweeps = radarDat.num_of_sweeps;
% update parameters
app2.bw = app2.fStop - app2.fStart;
app2.fc = app2.fStart + app2.bw/2;
app2.sampRate = 1/(app2.sTime)*app2.nSamp;

app2.prf = 1/(app2.sTime)/app2.nTx;
app2.rngMax = 3e8/(2*app2.bw)*app2.nSamp/2;
app2.velMax = app2.prf/2*3e8/app2.fc/2;

%%%%%%%%%%%%%%%%% read offline data into rawdata
%[filename,pathname] = uigetfile({'*.mat';'*.txt'});

%filename_p = strcat(pathname, filename);

% if strcmp(filename_p(end-3:end), '.mat')
%     radarDat = load(filename_p);
% 
%     app2.nRx = radarDat.nRx;
%     app2.nTx = radarDat.nTx;
%     app2.fStart = radarDat.fStart;
%     app2.fStop = radarDat.fStop;
%     app2.sTime = radarDat.sTime;%*1e-3;  
%     app2.nSamp = radarDat.nSamp;
%     rawdata = radarDat.ADCsamples;
% 
% elseif strcmp(filename_p(end-3:end), '.txt')
%     fileID = fopen(filename_p, 'r');
%     dataArray = textscan(fileID, '%f'); 
%     fclose(fileID);
% 
%     radarDat = dataArray{1};
%     app2.nRx = radarDat(1);
%     app2.nTx = radarDat(2);
%     app2.fStart = radarDat(3);
%     app2.fStop = radarDat(4);
%     app2.sTime = radarDat(5)*1e-3;
%     app2.nSamp = radarDat(6);
%     rawdata = radarDat(7:end);
% end

ind_Tx1 = find(rawdata>=49152);
ind_Tx1_diff = diff(ind_Tx1);
if any(ind_Tx1_diff ~= app2.nSamp*app2.nRx*2*app2.nTx) % *2: I&Q;
   fprintf("Some data were lost during USB data transmission.\n")
end

if app2.nTx == 2
    ind_Tx2 = find(rawdata>=32768 & rawdata < 49152);
    ind_Tx2_diff = diff(ind_Tx2);
    if any(ind_Tx2_diff ~= app2.nSamp*app2.nRx*2*app2.nTx) % *2: I&Q;
        fprintf("Some data were lost during USB data transmission.\n")
    end
end

if app2.nTx == 1
    rawdata(ind_Tx1) = rawdata(ind_Tx1) - 49152; % remove marker
    I1 = rawdata(ind_Tx1(1)+1:2*app2.nRx:ind_Tx1(end-1)+1+app2.nSamp); 
    Q1 = rawdata(ind_Tx1(1):2*app2.nRx:ind_Tx1(end-1)+app2.nSamp);
    complex_data = I1+1i*Q1;
%                     % Remove DC
%                     I1 = I1 - mean(I1);
%                     Q1 = Q1 - mean(Q1);
    if(app2.enable.HighpassEqui)
        complex_data = applyFilt(complex_data, app2);
    end 
    if app2.nRx >= 2
        I2 = rawdata(ind_Tx1(1)+3:2*app2.nRx:ind_Tx1(end-1)+3+app2.nSamp); 
        Q2 = rawdata(ind_Tx1(1)+2:2*app2.nRx:ind_Tx1(end-1)+2+app2.nSamp);
        complex_data = [complex_data, I2+1i*Q2];
%                         I2 = I2 - mean(I2);
%                         Q2 = Q2 - mean(Q2);
    end
    if app2.nRx == 4
        I3 = rawdata(ind_Tx1(1)+5:2*app2.nRx:ind_Tx1(end-1)+5+app2.nSamp); 
        Q3 = rawdata(ind_Tx1(1)+4:2*app2.nRx:ind_Tx1(end-1)+4+app2.nSamp); 
        I4 = rawdata(ind_Tx1(1)+7:2*app2.nRx:ind_Tx1(end-1)+7+app2.nSamp); 
        Q4 = rawdata(ind_Tx1(1)+6:2*app2.nRx:ind_Tx1(end-1)+6+app2.nSamp); 

        complex_data = [complex_data, I3+1i*Q3, I4+1i*Q4];
%                         I3 = I3 - mean(I3);
%                         Q3 = Q3 - mean(Q3);
%                         I4 = I4 - mean(I4);
%                         Q4 = Q4 - mean(Q4);
    end
end

if app2.nTx == 2
    ADCsamples_trimmed = rawdata(ind_Tx1(1):ind_Tx2(end-1)+app2.nSamp*2*app2.nRx-1);
    ADCsamples_m = reshape(ADCsamples_trimmed, app2.nSamp*2*app2.nRx, []);
    num_of_sweeps = size(ADCsamples_m,2);
    ADCsamples_t1_m = ADCsamples_m(:, 1:2:num_of_sweeps);
    ADCsamples_t2_m = ADCsamples_m(:, 2:2:num_of_sweeps);

    switch app2.nRx 
        case 4
            T1R1_m = 1i*ADCsamples_t1_m(1:8:app2.nSamp*8,:)+ADCsamples_t1_m(2:8:app2.nSamp*8+1,:);
            T1R1 = T1R1_m(:);
            T1R2_m = 1i*ADCsamples_t1_m(3:8:app2.nSamp*8+2,:)+ADCsamples_t1_m(4:8:app2.nSamp*8+3,:);
            T1R2 = T1R2_m(:);
            T1R3_m = 1i*ADCsamples_t1_m(5:8:app2.nSamp*8+4,:)+ADCsamples_t1_m(6:8:app2.nSamp*8+5,:);
            T1R3 = T1R3_m(:);
            T1R4_m = 1i*ADCsamples_t1_m(7:8:app2.nSamp*8+6,:)+ADCsamples_t1_m(8:8:app2.nSamp*8+7,:);
            T1R4 = T1R4_m(:);

            T2R1_m = 1i*ADCsamples_t2_m(1:8:app2.nSamp*8,:)+ADCsamples_t2_m(2:8:app2.nSamp*8+1,:);
            T2R1 = T2R1_m(:);
            T2R2_m = 1i*ADCsamples_t2_m(3:8:app2.nSamp*8+2,:)+ADCsamples_t2_m(4:8:app2.nSamp*8+3,:);
            T2R2 = T2R2_m(:);
            T2R3_m = 1i*ADCsamples_t2_m(5:8:app2.nSamp*8+4,:)+ADCsamples_t2_m(6:8:app2.nSamp*8+5,:);
            T2R3 = T2R3_m(:);
            T2R4_m = 1i*ADCsamples_t2_m(7:8:app2.nSamp*8+6,:)+ADCsamples_t2_m(8:8:app2.nSamp*8+7,:);
            T2R4 = T2R4_m(:);

            complex_data = [T1R1, T1R1];

        case 2

            T1R1_m = 1i*ADCsamples_t1_m(1:4:app2.nSamp*4,:)+ADCsamples_t1_m(2:4:app2.nSamp*4+1,:);
            T1R1 = T1R1_m(:);
            T1R2_m = 1i*ADCsamples_t1_m(3:4:app2.nSamp*4+2,:)+ADCsamples_t1_m(4:4:app2.nSamp*4+3,:);
            T1R2 = T1R2_m(:);

            T2R1_m = 1i*ADCsamples_t2_m(1:4:app2.nSamp*4,:)+ADCsamples_t2_m(2:4:app2.nSamp*4+1,:);
            T2R1 = T2R1_m(:);
            T2R2_m = 1i*ADCsamples_t2_m(3:4:app2.nSamp*4+2,:)+ADCsamples_t2_m(4:4:app2.nSamp*4+3,:);
            T2R2 = T2R2_m(:);
            
        %% Equiripple Highpass
        if(app2.enable.HighpassEqui)
            T1R1_f = applyFilt(T1R1, app2);
            T1R2_f = applyFilt(T1R2, app2);
            complex_data = [T1R1_f, T1R2_f];
        else
            complex_data = [T1R1, T1R2];
        end
            
        case 1
    end
end

end
function psd_plot(app,cmd_rd,X,Y,Z,clim_min,clim_max,xlab,ylab,ltx_fontSz)
%     imagesc(X, Y, Z);
    Z(Z <0.09) = 0;%eps;
    surf(X, Y, Z);view(2);
%     pcolor(X, Y, Z)%;view(2);
    %clim = get(gca,'CLim');
    %set(gca,'YDir','normal')
    %set(gca,'CLim',[clim_min, clim_max]);
    
    shading('interp')
    ylabel(ylab,'Interpreter','latex','Fontsize',ltx_fontSz)
    xlabel(xlab,'Interpreter','latex','Fontsize',ltx_fontSz)
    if strcmp(cmd_rd,'doppPSD')
        axis([0, X(end-1), -app.velMax, app.velMax])
    elseif strcmp(cmd_rd,'rdPSD')
        axis([0, X(end-1), -app.velMax, app.velMax])
    elseif strcmp(cmd_rd,'rngPSD')
        axis([0, X(end-1), Y(1), Y(end)])    
    end    
    c = colorbar; colormap((hot));%colormap(flipud(hot));
    c.TickLabelInterpreter = 'latex'; c.FontSize = ltx_fontSz;
    set(gca,'TickLabelInterpreter','latex','Fontsize',ltx_fontSz)
    drawnow  
%     exportgraphics(gca,'mysp600.pdf','Resolution',600)
end


function [dataVec_out] = applyFilt(dataVec, app2)


data_mat = reshape(dataVec(1:app2.nSamp*app2.nSweeps), app2.nSamp, []);
data_mat_f = HPFilter_2(data_mat,app2);

dataVec_out = data_mat_f(:);

end 

function [xf] = HPFilter_2(x,app2)

% ns = size(x,2);

% for i = 1:size(x,1)
%     xf(i,:)=filtfilt(signal_HPF,x(i, 1:ns));
% end

app2.signal_HPF = designfilt('highpassfir', 'StopbandFrequency',0.5, 'PassbandFrequency',3, 'StopbandAttenuation', 25, 'PassbandRipple', 1, 'SampleRate', app2.prf, 'DesignMethod', 'equiripple');
xf = transpose(filtfilt(app2.signal_HPF,transpose(x)));

% xf = transpose(diff(transpose(x)));

end