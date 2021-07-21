%% Description 
%  This script will observe a directory in a continuously and upon 
%  detecting a new file (recorded by the radar application), it will read 
%  the new file and process the data inside. 

%%  Remember to set the parameters for the radar the following
%-  Sweeptime = 1000
%-  Samples per sweep = 256
%-  Number of sweeps to process = 1000
%-  Must use the SISO configuration 
%-  These parameters are set by default. However, double check is
%-  important.

%% Requirements
%- Make sure if the cprintf toolbox is installed on the computer. This will
%- enable colored output.

%% Intrupt to exist from while loop.  
% ButtonHandle = uicontrol('Style', 'PushButton', ...
%                          'String', 'Stop Real-time Calssification Script', ...
%                          'Callback', 'delete(gcbf)', 'ForegroundColor', 'b', 'FontSize', 12);
% ButtonHandle.Position = [165 5 260 30]; 

%% 
clear all;
Nfftr = 512;
Nfftv = 512;

%% The path of the Folder to observe i.e., the location of the streamed data

folder_to_observe = strcat(pwd,'\ADC_streamedData\');

% Folder where spectral images will be saved. 
destDir1 = strcat(pwd,'\RealTime_Spectrograms\');

if (~exist(destDir1, 'dir'))
    mkdir(destDir1)
end 

dir_content = dir(folder_to_observe);
filenames = {dir_content.name};
current_files = filenames;

file_number = 0;

while true
    dir_content = dir(folder_to_observe);
    filenames = {dir_content.name};
    new_files = setdiff(filenames,current_files);
    if ~isempty(new_files)
        % deal with the new files
        current_files = filenames;
        pause(3)

        if(size(current_files,2)>2)
            
            fprintf('A new file is found\n')
            % New file will be the last file in the list
            filename = current_files{end};
            file_2_open = fullfile(folder_to_observe,filename);
            
            [fid, err] = fopen(file_2_open);
            
            if (fid >1)
                
                disp(filename)
                % Check how many links data is available    
                 [data, nSamp, sTime, nSweeps, fStart, fStop, nTx, nRx] =  getDataFromADC_Samples(file_2_open);
                
                % Check how many links data is available 
                
                fn = fieldnames(data);
                
                for k = 1:size(fn,1)                    
                    
                    %data_matrix = data.(fn{k});
                    [~] = computeStepsFromRangeDoppler(data, fStop, fStart, Nfftr, Nfftv, nTx, sTime, nSamp);
                    
                end 
                
                file_number = file_number+1;
                fclose('all');

            else
                
                pause(2)
                
                [fid, err] = fopen(file_2_open);
                
                % Check how many links data is available    
                 [data, nSamp, sTime, nSweeps, fStart, fStop, nTx, nRx] =  getDataFromADC_Samples(file_2_open);
                
                 % Check how many links data is available 
                 fn = fieldnames(data);
                              
                for k = 1:size(fn,1)
                    
                    %data_matrix = data.(fn{k});
                    
                    [~] = computeStepsFromRangeDoppler(data, fStop, fStart, Nfftr, Nfftv, nTx, sTime, nSamp);
                end
                file_number = file_number+1;
                fclose('all');
            end  
        end
    else
        
        fprintf('Waiting for a incomming data file...\n')
    end
    
%     if ~ishandle(ButtonHandle)
%         disp('User has stoped the Script');
%         
%         fclose('all');
%         break;
%     end
end
fclose('all');
