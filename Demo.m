
%% This script demonstrates the usage of RF step counter.

[file,path] = uigetfile;

%% Function call to get raw data from the RF data file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file = RF data file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Out Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data      = Data of an antenna link  
% fStop     = start of bandwidth at which data were recorder
% fStart    = end of bandwidth at which data were recorder
% nTx       = num of transmitting antennas 
% nRx       = num of receiving antennas 
% sTime     = Sweep Time   
% nSamp     = num of samples per sweep.
% nSweeps   = total number of sweeps

[data, nSamp, sTime, nSweeps, fStart, fStop, nTx, nRx] = getDataFromADC_Samples(fullfile(path, file));

%% Function call computeStepsFromRangeDoppler will perform following steps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 compute Range 
% 2 Compute Doppler Frequencies
% 3 Compute mean Doppler shift (MDS)
% 4 Count steps from MDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nfftr     = num of fft points for range
% Nfftv     = num of fft points for Doppler/velocity
% data      = Data of an antenna link  
% fStop     = start of bandwidth at which data were recorder
% fStart    = end of bandwidth at which data were recorder
% nTx       = num of transmitting antennas 
% nRx       = num of receiving antennas 
% sTime     = Sweep Time   
% nSamp     = num of samples per sweep.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Out param
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stepCount     = Number of steps in the file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
Nfftr = 1024; 
Nfftv = 1024;

[stepCount] = computeStepsFromRangeDoppler(data, fStop, fStart, Nfftr, Nfftv, nTx, sTime, nSamp);
toc;
fprintf('File %s contains %d steps\n', file, stepCount);
