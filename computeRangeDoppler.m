
function [stepCount] = computeRangeDoppler(data, fStop, fStart, Nfftr, Nfftv, nTx, sTime, nSamp)

%close all;
%clf;
data = data.T1R1;

data_matrix = (reshape(data, nSamp, []));

[bw, fc, sampRate, prf, rngMax, velMax, df_rng, df_dopp] = computeParam(fStop, fStart, Nfftr, Nfftv, nTx, sTime, nSamp);

%% Rym's Amplitude plots
% figure;
% t = linspace(0, size(data_matrix,2)/prf, size(data_matrix,2));
% plot(t, abs(sum(data_matrix,1)));
% xlabel('Time (s)', 'Interpreter','latex', 'FontSize',16);
% ylabel('Amplitude', 'Interpreter','latex', 'FontSize',16);
% drawnow;
%%

[Dopp_axis, Rng_axis] = computeRangeAxis (prf, fc, Nfftr, Nfftv, sTime, bw, sampRate);

[Rng_dft] = computeRangeProfiles(data_matrix, Nfftr, Rng_axis, prf);

[~, stepCount] = computeDopplerFromRange(Rng_dft, Nfftv, df_rng, df_dopp, prf, fc);

%computeVelocity(Rng_dft, Nfftv, velMax, df_rng, df_dopp, prf);

end 


%% This function computes parameters used in RangeProfile and DopplerRange, and velocity.
function [bw, fc, sampRate, prf, rngMax, velMax, df_rng, df_dopp] = computeParam(fStop, fStart, Nfftr, Nfftv, nTx, sTime, nSamp, nSweeps)

bw = fStop - fStart;

fc = fStart + bw/2;

sampRate = 1/(sTime)*nSamp;

prf= 1/(sTime)/nTx;

rngMax = 3e8/(2*bw)*nSamp/2;

%velMax = prf/2*3e8/fc/2;

velMax = 3e8/(4*fc*sTime);

df_rng = sampRate/Nfftr; % frequency resolution in range

df_dopp = prf/Nfftv;
end


%% This function computes Rangeprofiles
function [Rng_dft] = computeRangeProfiles(data_matrix, Nfftr, Rng_axis, prf)

data_matrix = bsxfun(@minus, data_matrix, mean(data_matrix, 2)); % remove static objects

% compute rangeProfiles 
Rng_dft = fft(data_matrix, Nfftr, 1);
Rng_dft = Rng_dft(1:Nfftr/2, :);



%%
%Rym's Amplitude plot
% figure;
% t = linspace(0, size(Rng_dft,2)/prf, size(Rng_dft,2));
% plot(t, abs(sum(Rng_dft,1)));
% xlabel('Time (s)', 'Interpreter','latex', 'FontSize',16);
% ylabel('Amplitude', 'Interpreter','latex', 'FontSize',16);
% drawnow;
%%


top = (abs(Rng_dft));
Rng_dft_norm = (top./abs(max(max(Rng_dft))));


% %plot Rangprofile
figure;
t = linspace(0, size(Rng_dft,2)/prf, size(Rng_dft,2)/prf);
imagesc(t, Rng_axis, 20*log10(Rng_dft_norm));
xlabel('Time (s)', 'Interpreter','latex', 'FontSize',16); 
ylabel('Range (m)', 'Interpreter','latex', 'FontSize',16); 


figure;

surf(20*log10(Rng_dft_norm (30:end, 8000:end-6000)));
shading interp;

end


%% This function computes RangeAxis
function [Dopp_axis, Rng_axis] = computeRangeAxis (prf, fc, Nfftr, Nfftv, sTime, bw, sampRate)

Dopp_axis = linspace(-prf/2,prf/2,Nfftv)*3e8/fc/2;
Rng_axis = linspace(0, sampRate/2, Nfftr/2+1)*3e8*sTime/(2*bw);

end

%% This function computes mean Doppler shift
function [mds] = computeMeanDopplerShift(Doppler, f)

top = sum(f'.*Doppler,1);
bottom = sum(Doppler, 1);
mds =  top./bottom;

end 

%% This function detects setps from the MDS
function [stepCount] = detectSteps(mds, t) 

if (mean(mds)<0)
    
    mds = mds*-1;
end

%figure;
hold on;
[pks,locs,~,~] = findpeaks(mds, t, 'threshold',0.001, 'MinPeakHeight',20, 'MinPeakDistance',0.005, 'MinPeakProminence',15);
plot(t,mds,locs, pks,'rv','MarkerFaceColor','r');
legend('Mean Doppler shift','Step marker');
count = num2str(size(pks,2));
stepTime = (mean(diff(locs)))*1000;
lgd = legend;
lgd.FontSize = 10;
lgd.Title.String = strcat('Total steps: ',count, string(newline), 'Mean step time (ms): ', num2str(stepTime));
grid 'on'
ax = gca;
ax.FontSize = 10;
xlabel('Time (s)', 'Interpreter','latex', 'FontSize',16);
ylabel('Frequency (Hz)', 'Interpreter','latex', 'FontSize',16);

stepCount = pks;
hold off;
end

%% This function computes Doppler from Range profiles.
function [RngDopp_psd, stepCount] = computeDopplerFromRange (Rng_dft, Nfftv, df_rng, df_dopp, prf, fc)

window = 64;
step = 1;
start = 1;
s_forward = 8;
Doppler = [];
tic;

while window <= size(Rng_dft, 2)
    
    windowed = conj(Rng_dft(:,start:window));
    
    windowed = bsxfun(@minus, windowed, mean(windowed, 2)); % remove stationary objects
    
    %jj = abs(fft(windowed,Nfftv,2)).^2;
    %kk = abs(fftshift(jj(1:Nfftv/2, :),2)).^2;
    
    RngDopp_psd = abs(fftshift(fft(windowed,Nfftv,2),2)).^2;%/(df_dopp*df_rng); %Works fine without this division
    
    Doppler(:,step) = (sum(RngDopp_psd,1))/max(sum(RngDopp_psd,1));
    
    step = step+1;
    
    start = start+s_forward;
    window = window+s_forward;
    
end

% Plot Doppler w.r.t time

figure;

vel_scale = linspace(-prf/2,prf/2,Nfftv)*3e8/fc/2;
dopp_scale  =  dopplerFromSpeed(fc, vel_scale);

t = linspace(0, size(Rng_dft,2)/prf, size(Doppler,2));
imagesc(t, dopp_scale, flipud(Doppler));
set(gca,'YDir','normal')
%ylim([-Nfftv/2 Nfftv/2])
view(2);
shading interp;
xlabel('Time (s)', 'Interpreter','latex', 'FontSize',16);
ylabel('Frequency (Hz)', 'Interpreter','latex', 'FontSize',16);
drawnow;
toc;

%% Compute mean Doppler shift.

[mds] = computeMeanDopplerShift(Doppler, -dopp_scale);

mds_f =  sgolayfilt(mds, 3, 51);
%mds((mds < 20) & (mds > -20)) =0;

figure;
plot(t, mds_f);

%[stepCount] = detectSteps(mds_f, t); 


%speedFromDoppler_2(Doppler, -vel_scale, t);

end


function[] = speedFromDoppler_2(Doppler, vel_scale, t)

top = sum(vel_scale'.*Doppler,1);
bottom = sum(Doppler, 1);
vel =  top./bottom;

figure;
plot(t, vel);

xlabel('Time (s)', 'Interpreter','latex', 'FontSize',16); 
ylabel('Velocity (m/s)', 'Interpreter','latex', 'FontSize',16); 

end 

%% This function computes velocity from the range profile.

function computeVelocity(Rng_dft, Nfftv, velMax, df_rng, df_dopp, prf)

window = 32;
step = 1;
start = 1;
s_forward = 8;

vel_axis = linspace(-velMax, velMax, Nfftv);

while window <= size(Rng_dft, 2)

windowed = Rng_dft(:,start:window);
    
RngDopp_psd = abs(fftshift(fft(windowed,Nfftv,2),2)).^2/(df_dopp*df_rng); %Works fine without this division

[Doppler_psd_max, Doppler_psd_max_index] = max(max(RngDopp_psd));

vel(:,step) = vel_axis(Doppler_psd_max_index);

step = step+1;
       
start = start+s_forward;
window = window+s_forward;

end 

% Velocity plot
figure;
t = linspace(0, size(Rng_dft,2)/prf, length(vel));
plot(t, sgolayfilt(vel, 3, 11));
xlabel('Time (s)', 'Interpreter','latex', 'FontSize',16); 
ylabel('Velocity (m/s)', 'Interpreter','latex', 'FontSize',16); 

end 



function [untrended] =  untrend(t, mds)
opol = 2;
[p,s,mu] = polyfit(t,mds,opol);
f_y = polyval(p,t,[],mu);

untrended = mds - f_y;

end 


function [xf] = HPFilter(x)

 ns = size(x,2);
[b,a] = butter(3,0.005, 'high');
[h, f1] = freqz(b, a, ns);

xf = filtfilt(b,a,x(1:ns));

end 

function [pcComponents] =  getMinorPCAcomponents(x)

[coeff, score, ~, ~, explained, ~] = pca(x');

pcComponents = score';

end 

function vel =  speedFromDoppler(fc, dopshift, t)

lambda = physconst('LightSpeed')/fc/2;

vel = dop2speed(sgolayfilt(dopshift, 3, 11),lambda); 
figure;
plot(t, vel);

xlabel('Time (s)', 'Interpreter','latex', 'FontSize',16); 
ylabel('Velocity (m/s)', 'Interpreter','latex', 'FontSize',16); 

end 

function Dopp =  dopplerFromSpeed(fc, vel)

lambda = physconst('LightSpeed')/fc/2;
Dopp = speed2dop(vel,lambda); 

end 



