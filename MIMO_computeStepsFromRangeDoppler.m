%% This function will first compute Doppler from Rnage profiles. 

%% Inputs of all functions:

% data      = Data of an antenna link  
% fStop     = start of bandwidth at which data were recorder
% fStart    = end of bandwidth at which data were recorder
% Nfftr     = # of fft points for range
% Nfftv     = # of fft points for Doppler/velocity
% nTx       = # of transmitting antennas 
% sTime     = Sweep Time   
% nSamp     = # of samples per sweep.
% data_matrix = Reshaped time series data to a matrix
% Rng_axis  = Range axis
% prf       = pulse resolution frequency 
% bw        = bandwidth
% sampRate  = sampling rate 
% rngMax    = Maximum range
% velMax    = maximum velocity
% df_rng    = frequency resolution in range
% df_dopp   = frequency resolution in doppler

function [stepCount] = computeStepsFromRangeDoppler(data, fStop, fStart, Nfftr, Nfftv, nTx, sTime, nSamp)
close all;

% data      = Data of an antenna link  
% fStop     = start of bandwidth at which data were recorder
% fStart    = end of bandwidth at which data were recorder
% Nfftr     = # of fft points for range
% Nfftv     = # of fft points for Doppler/velocity
% nTx       = # of transmitting antennas 
% sTime     = Sweep Time   
% nSamp     = # of samples per sweep.


% Reshape time series data to a matrix
data = data.T1R2;
data_matrix = (reshape(data, nSamp, []));

% First compute some basic parameters 
[bw, fc, sampRate, prf, rngMax, velMax, df_rng, df_dopp] = computeParam(fStop, fStart, Nfftr, Nfftv, nTx, sTime, nSamp);

% Compute Range and Doppler axes
[Dopp_axis, Vel_axis, Rng_axis] = computeRangeAxis (prf, fc, Nfftr, Nfftv, sTime, bw, sampRate);

[Rng_dft] = computeRangeProfiles(data_matrix, Nfftr, Rng_axis, prf);

[Doppler, mds, t] =  computeDopplerFromRange(Rng_dft, Nfftv, df_rng, df_dopp, prf, fc, Dopp_axis, Rng_axis);

[stepCount] = detectSteps(mds, t);
stepCount =0;

%% The following two methods are for deomastration that average or summing delays to compute CTF and then 
 % use CTF to compute Doppler introduce noise when there is no activity. Consequently causing huge 
 % fluctuations in the mean Doppler shift.
 
%[STFT, mds_f, t] = computeDopplerFromRange_CTF (Rng_dft, Nfftv, df_rng, df_dopp, prf, fc, Dopp_axis);
%[Doppler1, mds_f, t] = computeDopplerFromRange_CTF_FFT (Rng_dft, Nfftv, df_rng, df_dopp, prf, fc, Dopp_axis);

end 

%% This function computes parameters used in RangeProfile and DopplerRange, and velocity.
function [bw, fc, sampRate, prf, rngMax, velMax, df_rng, df_dopp] = computeParam(fStop, fStart, Nfftr, Nfftv, nTx, sTime, nSamp, nSweeps)

% bandwidth
bw = fStop - fStart;

% Central carrier frequency
fc = fStart + bw/2;

% Sampling rate
sampRate = 1/(sTime)*nSamp;

% Pulse rsoulution freq
prf= 1/(sTime)/nTx;

% Maximum range
rngMax = physconst('LightSpeed')/(2*bw)*nSamp/2;

% Maximum velocity
%velMax = prf/2*3e8/fc/2;
velMax = physconst('LightSpeed')/(4*fc*sTime);

% frequency resolution in range
df_rng = sampRate/Nfftr; 

% frequency resolution in Doppler
df_dopp = prf/Nfftv;
end


%% This function computes RangeAxis
function [Dopp_axis, Vel_axis, Rng_axis] = computeRangeAxis (prf, fc, Nfftr, Nfftv, sTime, bw, sampRate)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prf       = pulse resolution frequency 
% Nfftr     = # of fft points for range
% Nfftv     = # of fft points for Doppler/velocity
% sTime     = Sweep Time  
% bw        = bandwidth
% sampRate  = sampling rate 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rng_axis  = Range axis
% Dopp_axis = Doppler axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vel_axis = linspace(-prf/2,prf/2,Nfftv)*3e8/fc/2;
Rng_axis = linspace(0, sampRate/2, Nfftr/2+1)*3e8*sTime/(2*bw);
Dopp_axis =  dopplerFromSpeed(fc, Vel_axis);
end

%% This function computes Rangeprofiles
function [Rng_dft] = computeRangeProfiles(data_matrix, Nfftr, Rng_axis, prf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data_matrix = Reshaped time series data to a matrix
% Nfftr     = # of fft points for range
% Rng_axis  = Range axis
% prf       = pulse resolution frequency 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rng_dft   = Range Profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove static objects

[data_matrix1] = HPFilter_2(data_matrix);

data_matrix1 = bsxfun(@minus, data_matrix1, mean(data_matrix1, 2)); 

% Compute rangeProfiles 
Rng_dft = fft(data_matrix1, Nfftr, 1);
Rng_dft = Rng_dft(1:Nfftr/2, :);

% Normalize Range for plotting
top = abs(Rng_dft).^2;
Rng_dft_norm = top./sum(top,1);

%plot Rangprofile
plotting = 1;
if(plotting ==1)
    figure;
    t = linspace(0, size(Rng_dft,2)/prf, size(Rng_dft,2)/prf);
    colormap(hot);
    imagesc(t, Rng_axis, Rng_dft_norm);
    xlabel('Time, t (s)', 'Interpreter','latex', 'FontSize',16);
    ylabel('Range (m)', 'Interpreter','latex', 'FontSize',16);
end


% Restricted range plot
%range = 35;

%[i, j, ~]=  find(Rng_axis <= range);

% if(plotting ==1)
%     figure;
%     t = linspace(0, size(Rng_dft,2)/prf, size(Rng_dft,2)/prf);
%     colormap(hot);
%     imagesc(t, Rng_axis(j(1):j(end)), Rng_dft_norm(j(1):j(end), :));
%     xlabel('Time, t (s)', 'Interpreter','latex', 'FontSize',16);
%     ylabel('Range (m)', 'Interpreter','latex', 'FontSize',16);
% end

end


%% Compute CTF from range and then compute spectrogram.
function [norm_spectrogram, mds_f, t] = computeDopplerFromRange_CTF (Rng_dft, Nfftv, df_rng, df_dopp, prf, fc, Dopp_axis)
Fs = prf;

Rng_dft = bsxfun(@minus, Rng_dft, mean(Rng_dft, 2)); 

tf = conj(sum(Rng_dft,1));

[STFT,f,t] = spectrogram(tf,rectwin(64),56,Nfftv,Fs,'centered','yaxis');

Sxx1 = abs(STFT).^2;

norm_spectrogram = Sxx1./sum(Sxx1, 1);

norm_spectrogram(norm_spectrogram <0.005) = 0;%eps;
figure;
colormap(hot)
%t = linspace(0, size(Rng_dft,2)/prf, size(Doppler,2));
imagesc(t, Dopp_axis, norm_spectrogram);
set(gca,'YDir','normal')
%ylim([-Nfftv/2 Nfftv/2])
view(2);
shading interp;
xlabel('Time (s)', 'Interpreter','latex', 'FontSize',16);
ylabel('Frequency (Hz)', 'Interpreter','latex', 'FontSize',16);
drawnow;

%% Compute mean Doppler shift.
f = f'; % so that next function works correctly
[mds] = computeMeanDopplerShift(norm_spectrogram, f);
mds_f =  sgolayfilt(mds, 3, 51);
figure;
plot(t, mds_f);
xlabel('Time, t (s)', 'Interpreter','latex', 'FontSize',16);
ylabel('Frequency f, (Hz)', 'Interpreter','latex', 'FontSize',16);
end 
%% This function computes Doppler from Range profiles.
function [Doppler, mds_f, t] = computeDopplerFromRange(Rng_dft, Nfftv, df_rng, df_dopp, prf, fc, Dopp_axis, Rng_axis)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rng_dft   = Range Profile
% Nfftv     = # of fft points for Doppler/velocity
% df_rng    = frequency resolution in range
% df_dopp   = frequency resolution in Doppler
% prf       = pulse resolution frequency 
% fc        = Central carrier frequency
% Dopp_axis = Doppler Axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dopp_psd  = Doppler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Rng_dft = bsxfun(@minus, Rng_dft, mean(Rng_dft, 2)); 


%% Range sets 
R11_start = max(find(Rng_axis<=0.6));
R11_end = max(find(Rng_axis<=5.8));


RInt_start = max(find(Rng_axis<=5.85));
RInt_end = max(find(Rng_axis<=11.66));


R22_start = max(find(Rng_axis<=12.0));
R22_end = max(find(Rng_axis<=20.00));


ranges= [R11_start, R11_end;RInt_start,RInt_end; R22_start,R22_end];


%%

for i = 1:size(ranges,1)
window = 64;
step = 1;
start = 1;
s_forward = 8;
Doppler = [];
    while window <= size(Rng_dft, 2)
        
        windowed = conj(Rng_dft(ranges(i,1):ranges(i,2),start:window));
        % remove stationary objects
        windowed = bsxfun(@minus, windowed, mean(windowed, 2));
        
        RngDopp_psd = abs(fftshift(fft(windowed,Nfftv,2),2)).^2/(df_dopp*df_rng);
        Doppler(:,step) = sum(RngDopp_psd,1)./sum(sum(RngDopp_psd,1));
        
        step = step+1;
        start = start+s_forward;
        window = window+s_forward;
        
    end
    window = 0;
    % Threshold the Doppler
    Doppler(Doppler <0.005) = 0;%eps;
    %Plot Doppler w.r.t time
    
    t = linspace(0, size(Rng_dft,2)/prf, size(Doppler,2));
    plotting = 1;
    if(plotting ==1)
        figure;
        colormap(hot);colorbar;
        surf(t, Dopp_axis,((Doppler)));
        set(gca,'YDir','normal')
        %ylim([-Nfftv/2 Nfftv/2])
        view(2);
        shading interp;
        xlabel('Time, t (s)', 'Interpreter','latex', 'FontSize',16);
        ylabel('Doppler frequency, f (Hz)', 'Interpreter','latex', 'FontSize',16);
        drawnow;
    end
    
    %% Compute mean Doppler shift.
    [mds] = computeMeanDopplerShift(Doppler, Dopp_axis);
    mds(isnan(mds))=0;
    mds_f =  sgolayfilt(mds, 3, 71);
    plotting = 1;
    if(plotting ==1)
        figure;
        plot(t, mds_f);
        xlabel('Time, t (s)', 'Interpreter','latex', 'FontSize',16);
        ylabel('Mean Dopper shift, $B_f(t)$ (Hz)', 'Interpreter','latex', 'FontSize',16);
        grid on;
        drawnow;
    end
end

end

%% Compute CTF from range and then compute 
function [Doppler, mds_f, t] = computeDopplerFromRange_CTF_FFT (Rng_dft, Nfftv, df_rng, df_dopp, prf, fc, Dopp_axis)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rng_dft   = Range Profile
% Nfftv     = # of fft points for Doppler/velocity
% df_rng    = frequency resolution in range
% df_dopp   = frequency resolution in Doppler
% prf       = pulse resolution frequency 
% fc        = Central carrier frequency
% Dopp_axis = Doppler Axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dopp_psd  = Doppler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

window = 64;
step = 1;
start = 1;
s_forward = 8;
Doppler = [];



tf = conj(sum(Rng_dft,1));

while window <= size(tf, 2)
    
    windowed = (tf(:,start:window));
    % remove stationary objects
    windowed = bsxfun(@minus, windowed, mean(windowed, 2)); 
    
    RngDopp_psd = abs(fftshift(fft(windowed,Nfftv,2),2)).^2/(df_dopp*df_rng);
    Doppler(:,step) = sum(RngDopp_psd,1)./sum(sum(RngDopp_psd,1));
    
    step = step+1;
    start = start+s_forward;
    window = window+s_forward;
    
end

% Threshold the Doppler
Doppler(Doppler <0.095) = 0;%eps;
%Plot Doppler w.r.t time
figure;
colormap(hot)
t = linspace(0, size(Rng_dft,2)/prf, size(Doppler,2));
imagesc(t, Dopp_axis, Doppler);
set(gca,'YDir','normal')
%ylim([-Nfftv/2 Nfftv/2])
view(2);
shading interp;
xlabel('Time, t (s)', 'Interpreter','latex', 'FontSize',16);
ylabel('Doppler frequency, f (Hz)', 'Interpreter','latex', 'FontSize',16);
drawnow;

%% Compute mean Doppler shift.
[mds] = computeMeanDopplerShift(Doppler, Dopp_axis);
mds_f =  sgolayfilt(mds, 3, 51);
figure;
plot(t, mds_f);
xlabel('Time, t (s)', 'Interpreter','latex', 'FontSize',16);
ylabel('Mean Dopper shift, $B_f(t)$ (Hz)', 'Interpreter','latex', 'FontSize',16);
end

%% This function computes mean Doppler shift
function [mds] = computeMeanDopplerShift(Doppler, f)
mds =  sum(f'.*Doppler,1)./sum(Doppler, 1);
end 

%% Function to convert velocity to Doppler
function Dopp =  dopplerFromSpeed(fc, vel)
lambda = physconst('LightSpeed')/fc/2;
Dopp = speed2dop(vel,lambda); 
end 

%% This function detects setps from the MDS
function [stepCount] = detectSteps(mds, t) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mds    = mean Doppler Shift
%t      = time 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stepCount = number of steps reurned
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert MDS to only positive side so that we only need to detect peaks.
if (mean(mds)<0)
    mds = mds*-1;
end

[pks,locs,~,~] = findpeaks(mds, t);%findpeaks(mds, t, 'threshold',0.001, 'MinPeakHeight',20, 'MinPeakDistance',0.005, 'MinPeakProminence',15);
plotting = 1;
if(plotting == 1)
    figure;
    hold on;
    plot(t,mds,locs, pks,'rv','MarkerFaceColor','r');
    legend('Mean Doppler shift','Step marker');
    count = num2str(size(pks,2));
    stepTime = (mean(diff(locs)))*1000;
    lgd = legend;
    lgd.FontSize = 10;
    lgd.Title.String = strcat('Total steps: ',count, string(newline));%, 'Mean step time (ms): ', num2str(stepTime));
    grid 'on'
    ax = gca;
    ax.FontSize = 10;
    xlabel('Time, t (s)', 'Interpreter','latex', 'FontSize',16);
    ylabel('Mean Dopper shift, $B_f(t) (Hz)$', 'Interpreter','latex', 'FontSize',16);
end
stepCount = size(pks,2);
hold off;
end

function [xf] = HPFilter(x)

ns = size(x,2);
[b,a] = butter(3,0.005, 'high');
[h, f1] = freqz(b, a, ns);

for i = 1:size(x,1)
    
    xf(i,:)=filtfilt(b,a,x(i, 1:ns));
    
end

end

function [xf] = HPFilter_2(x)

ns = size(x,2);

signal_HPF = designfilt('highpassfir', 'StopbandFrequency',0.5, 'PassbandFrequency',15, 'StopbandAttenuation', 15, 'PassbandRipple', 1, 'SampleRate', 1000, 'DesignMethod', 'equiripple');

for i = 1:size(x,1)
    xf(i,:)=filtfilt(signal_HPF,x(i, 1:ns));
end

end



