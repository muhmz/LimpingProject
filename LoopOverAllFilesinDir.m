clear all;
close all;

dataDir = uigetdir; %gets directory

dataFiles = dir(fullfile(dataDir,'*.mat')); %gets all mat files in struct

data_output_dir = strcat(pwd,'/Data/Output/6CM_Right_leg/');

mkdir(strcat(pwd,'/Data/Output/6CM_Right_leg/'));

file_name = {};
steps = [];
ds ={};
locs = {};
for k = 1:length(dataFiles)
    
    baseFileName = dataFiles(k).name;
    
    fullFileName = fullfile(dataDir, baseFileName);
    
    fprintf(1, 'Now processing %s\n', fullFileName);
    
    
    %Call your processing code
    [data, nSamp, sTime, nSweeps, fStart, fStop, nTx, nRx] = getDataFromADC_Samples(fullFileName);
    
    Nfftr = 1024;
    Nfftv = 1024;
    
    
    [stepCount, mds, peaks_locs] = computeStepsFromRangeDoppler(data, fStop, fStart, Nfftr, Nfftv, nTx, sTime, nSamp);
    
    file_name(k,:) = cellstr(baseFileName);
    
    steps = [steps; stepCount];
    
    ds(k,:) = {mds};
    
    locs(k,:) = {peaks_locs};
    
end

fileName = file_name;
results.fileName = fileName;
results.steps = steps;
results.ds = ds;
results.locs = locs;

save(strcat(data_output_dir,'\results_S1_6CM_Right.mat'), 'results');
