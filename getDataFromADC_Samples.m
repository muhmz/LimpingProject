function [data, nSamp, sTime, nSweeps, fStart, fStop, nTx, nRx] = getDataFromADC_Samples(file)

% This function will read data from the mat file 

radarDat = load(file);

nRx = radarDat.nRx;
nTx = radarDat.nTx;
fStart = radarDat.fStart;
fStop = radarDat.fStop;
sTime = radarDat.sTime;
nSamp = radarDat.nSamp;
%nSweeps = radarDat.num_of_sweeps;
ADCsamples = radarDat.ADCsamples;

% Now convert ADC Samples to complex data for each link

%%
% A larger number, '49152', was added to the first sample from Q channel to identify 
% the beginning of each sweep, so that we can check if there is any data
% loss in each sweep

ind_Tx1 = find(ADCsamples>=49152);
ind_Tx1_diff = diff(ind_Tx1);
if any(ind_Tx1_diff ~= nSamp*nRx*2*nTx) % *2: I&Q;
   fprintf("Some data were lost during USB data transmission.\n")
end

if nTx == 2
    ind_Tx2 = find(ADCsamples>=32768 & ADCsamples < 49152);
    ind_Tx2_diff = diff(ind_Tx2);
    if any(ind_Tx2_diff ~= nSamp*nRx*2*nTx) % *2: I&Q;
        fprintf("Some data were lost during USB data transmission.\n")
    end
end

%%
fc = (fStart + fStop)/2; % center frequency
Bandwidth = fStop - fStart; % in Hz
Samp_rate = nSamp/sTime;
T =  length(ADCsamples)/nRx/2/Samp_rate; % total observation time
nSweeps = T/(sTime); % number of sweeps
F = 1/sTime; % Doppler span -F/2 -> F/2

%% Truncate additional sweeps


%%
% The recorded ADC samples of each sweep are formatted as follow:
%
% 1) For 1 Tx 1 Rx models
% I(1), Q(1) + 49152, I(2), Q(2), I(3), Q(3), ..., I(nSamp), Q(nSamp)
% Q(2), I(2), Q(3), I(3), ..., Q(nSamp), I(nSamp)
% ... (nSweeps)
%
% 2) For 1 Tx 2 Rx models
% I1(1), Q1(1) + 49152, I2(1), Q2(1), I1(2), Q1(2), I2(2), Q2(2), ...,
% I1(nSamp), Q1(nSamp), I2(nSamp), Q2(nSamp)
% ... (nSweeps)
%
% 3) For 1 Tx 4 Rx models
% I1(1), Q1(1) + 49152, I2(1), Q2(1), I3(1), Q3(1), I4(1), Q4(1),
% I1(2), Q1(2), I2(2), Q2(2), I3(2), Q3(2), I4(2), Q4(2),
% ...
% I1(nSamp), Q1(nSamp), I2(nSamp), Q2(nSamp), I3(nSamp), Q3(nSamp), I4(nSamp), Q4(nSamp)
% ... (nSweeps)

% In extracting the I and Q data, we need to hunt for the first marker and rearrange the ADC samples.
if nTx == 1
     marker_1st = ind_Tx1(1);
     marker_last = ind_Tx1(end-1);
%     I1 = ADCsamples(marker_1st+1:2*nRx:marker_last+1+nSamp); 
%     Q1 = ADCsamples(marker_1st:2*nRx:marker_last+nSamp); 
    ADCsamples(ind_Tx1) = ADCsamples(ind_Tx1) - 49152; % remove marker
    I1 = ADCsamples(ind_Tx1(1)+1:2*nRx:ind_Tx1(end-1)+1+nSamp); 
    Q1 = ADCsamples(ind_Tx1(1):2*nRx:ind_Tx1(end-1)+nSamp);
    % Remove DC
    I1 = I1 - mean(I1);
    Q1 = Q1 - mean(Q1);
    if nRx >= 2
        I2 = ADCsamples(ind_Tx1(1)+3:2*nRx:ind_Tx1(end-1)+3+nSamp); 
        Q2 = ADCsamples(ind_Tx1(1)+2:2*nRx:ind_Tx1(end-1)+2+nSamp);
        I2 = I2 - mean(I2);
        Q2 = Q2 - mean(Q2);
    end
    if nRx == 4
        I3 = ADCsamples(ind_Tx1(1)+5:2*nRx:ind_Tx1(end-1)+5+nSamp); 
        Q3 = ADCsamples(ind_Tx1(1)+4:2*nRx:ind_Tx1(end-1)+4+nSamp); 
        I4 = ADCsamples(ind_Tx1(1)+7:2*nRx:ind_Tx1(end-1)+7+nSamp); 
        Q4 = ADCsamples(ind_Tx1(1)+6:2*nRx:ind_Tx1(end-1)+6+nSamp); 
        I3 = I3 - mean(I3);
        Q3 = Q3 - mean(Q3);
        I4 = I4 - mean(I4);
        Q4 = Q4 - mean(Q4);
    end
end

% for 2Tx 4Rx models
% a) in 2Tx 4Rx mode:
% I_t1r1(1), Q_t1r1(1)+46152, I_t1r2(1), Q_t1r2(1), I_t1r3(1), Q_t1r3(1),I_t1r4(1), Q_t1r4(1),
% I_t1r1(2), Q_t1r1(2), I_t1r2(2), Q_t1r2(2), I_t1r3(2), Q_t1r3(2),I_t1r4(2), Q_t1r4(2)
% ...
% I_t1r1(nSamp), Q_t1r1((nSamp), I_t1r2(nSamp), Q_t1r2(nSamp), I_t1r3(nSamp), Q_t1r3(nSamp), I_t1r4(nSamp), Q_t1r4(nSamp),
% I_t2r1(1), Q_t2r1(1)+32768, I_t2r2(1), Q_t2r2(1), I_t2r3(1), Q_t2r3(1),I_t2r4(1), Q_t2r4(1)
% I_t2r1(2), Q_t2r1(2), I_t2r2(2), Q_t2r2(2), I_t2r3(2),Q_t2r3(2),I_t2r4(2), Q_t2r4(2),
% ...
% I_t2r1(nSamp), Q_t2r1((nSamp), I_t2r2(nSamp), Q_t2r2(nSamp), I_t2r3(nSamp), Q_t2r3(nSamp), I_t2r4(nSamp), Q_t2r4(nSamp),
% ... (nSweeps)

% b) in 2Tx 2Rx mode:
% I_t1r1(1), Q_t1r1(1)+46152, I_t1r2(1), Q_t1r2(1),  
% I_t1r1(2), Q_t1r1(2), I_t1r2(2), Q_t1r2(2), 
% ...
% I_t1r1(nSamp), Q_t1r1((nSamp), I_t1r2(nSamp), Q_t1r2(nSamp),
% I_t2r1(1), Q_t2r1(1)+32768, I_t2r2(1), Q_t2r2(1), 
% I_t2r1(2), Q_t2r1(2), I_t2r2(2), Q_t2r2(2),
% ...
% I_t2r1(nSamp), Q_t2r1((nSamp), I_t2r2(nSamp), Q_t2r2(nSamp), 
% ... (nSweeps)

% c) in 2Tx 1Rx mode:
% I_t1r1(1), Q_t1r1(1)+46152,
% I_t1r1(2), Q_t1r1(2),  
% ...
% I_t1r1(nSamp), Q_t1r1((nSamp), 
% I_t2r1(1), Q_t2r1(1)+32768, 
% I_t2r1(2), Q_t2r1(2), 
% ...
% I_t2r1(nSamp), Q_t2r1((nSamp), 
% ... (nSweeps)
marker_1st = ind_Tx1(1);

if nTx == 2
    marker_last = ind_Tx2(end);
    ADCsamples_trimmed = ADCsamples(ind_Tx1(1):ind_Tx2(end-1)+nSamp*2*nRx-1);
    ADCsamples_m = reshape(ADCsamples_trimmed, nSamp*2*nRx, []);
    num_of_sweeps = size(ADCsamples_m,2);
    ADCsamples_t1_m = ADCsamples_m(:, 1:2:num_of_sweeps);
    ADCsamples_t2_m = ADCsamples_m(:, 2:2:num_of_sweeps);
    
    switch nRx 
        case 4
            T1R1_m = 1i*ADCsamples_t1_m(1:8:nSamp*8,:)+ADCsamples_t1_m(2:8:nSamp*8+1,:);
            T1R1 = T1R1_m(:);%,1:floor(nSweeps/2));
            T1R2_m = 1i*ADCsamples_t1_m(3:8:nSamp*8+2,:)+ADCsamples_t1_m(4:8:nSamp*8+3,:);
            T1R2 = T1R2_m(:);%,1:floor(nSweeps/2));
            T1R3_m = 1i*ADCsamples_t1_m(5:8:nSamp*8+4,:)+ADCsamples_t1_m(6:8:nSamp*8+5,:);
            T1R3 = T1R3_m(:);%,1:floor(nSweeps/2));
            T1R4_m = 1i*ADCsamples_t1_m(7:8:nSamp*8+6,:)+ADCsamples_t1_m(8:8:nSamp*8+7,:);
            T1R4 = T1R4_m(:);%,1:floor(nSweeps/2));
            
            T2R1_m = 1i*ADCsamples_t2_m(1:8:nSamp*8,:)+ADCsamples_t2_m(2:8:nSamp*8+1,:);
            T2R1 = T2R1_m(:);%,1:floor(nSweeps/2));
            T2R2_m = 1i*ADCsamples_t2_m(3:8:nSamp*8+2,:)+ADCsamples_t2_m(4:8:nSamp*8+3,:);
            T2R2 = T2R2_m(:);%,1:floor(nSweeps/2));
            T2R3_m = 1i*ADCsamples_t2_m(5:8:nSamp*8+4,:)+ADCsamples_t2_m(6:8:nSamp*8+5,:);
            T2R3 = T2R3_m(:);%,1:floor(nSweeps/2));
            T2R4_m = 1i*ADCsamples_t2_m(7:8:nSamp*8+6,:)+ADCsamples_t2_m(8:8:nSamp*8+7,:);
            T2R4 = T2R4_m(:);%,1:floor(nSweeps/2));
            
        case 2
            
            T1R1_m = 1i*ADCsamples_t1_m(1:4:nSamp*4,:)+ADCsamples_t1_m(2:4:nSamp*4+1,:);
            T1R1 = T1R1_m(:);%,1:floor(nSweeps/2));
            T1R2_m = 1i*ADCsamples_t1_m(3:4:nSamp*4+2,:)+ADCsamples_t1_m(4:4:nSamp*4+3,:);
            T1R2 = T1R2_m(:);%,1:floor(nSweeps/2));
            
            T2R1_m = 1i*ADCsamples_t2_m(1:4:nSamp*4,:)+ADCsamples_t2_m(2:4:nSamp*4+1,:);
            T2R1 = T2R1_m(:);%,floor(1:nSweeps/2));
            T2R2_m = 1i*ADCsamples_t2_m(3:4:nSamp*4+2,:)+ADCsamples_t2_m(4:4:nSamp*4+3,:);
            T2R2 = T2R2_m(:);%,floor(1:nSweeps/2));
                     
            
        case 1
            
    end
end


%% preparing return structure...

if (nTx == 1 && nRx == 1)
    
    complexData_1 = complex(I1 ,Q1);
    
    complexData_1 = complexData_1(1:floor(nSamp*(nSweeps-4)), :);
    data = struct('T1R1', complexData_1);
    
end 

if (nTx == 1 && nRx == 2)
    
    complexData_1 = complex(I1 , Q1);
    complexData_2 = complex(I2 , Q2);
    
    
    complexData_1 = complexData_1(1:floor(nSamp*(nSweeps-4)), :);
    complexData_2 = complexData_2(1:floor(nSamp*(nSweeps-4)), :);

    
    data = struct('T1R1', complexData_1, 'T1R2', complexData_2);
        
end 

if (nTx == 1 && nRx == 4)
    
    complexData_1 = complex(I1 , Q1);
    complexData_2 = complex(I2 , Q2);
    complexData_3 = complex(I3 , Q3);
    complexData_4 = complex(I4 , Q4);
    
    
    
    complexData_1 = complexData_1(1:floor(nSamp*(nSweeps-4)), :);
    complexData_2 = complexData_2(1:floor(nSamp*(nSweeps-4)), :);
    complexData_3 = complexData_3(1:floor(nSamp*(nSweeps-4)), :);
    complexData_4 = complexData_4(1:floor(nSamp*(nSweeps-4)), :);
    
    data = struct('T1R1', complexData_1, 'T1R2', complexData_2, 'T1R3',complexData_3 , 'T1R4', complexData_4);
    
end

if(nTx ==2 && nRx ==4)
    
    %data_combined= [T1R1;T1R2;T1R3;T1R4;T2R1;T2R2;T2R3;T2R4];
    %data = struct('T_All', data_combined);
    
    data = struct('T1R1', T1R1, 'T1R2', T1R2, 'T1R3', T1R3 , 'T1R4', T1R4, 'T2R1', T2R1, 'T2R2', T2R2, 'T2R3', T2R3 , 'T2R4', T2R4);

end


if(nTx ==2 && nRx ==2)
    
    data = struct('T1R1', T1R1, 'T1R2', T1R2, 'T2R1', T2R1, 'T2R2', T2R2);

end

end 

