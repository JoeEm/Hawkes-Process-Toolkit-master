% function [Seqs, Stats] = LoadData(file_path, file_format, ...
%                                 time_format, time_offset, time_scale)
%                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Load data and convert it to matlab format
%
% file_path: the path of data
% file_format: .csv or .txt file
% time_format: currently we support two formats. 1) real number; 2) real
% time, e.g., "year/month/day hour:minute:second" 3) ... you can define
% your format
%
%
% Remember to reverse the List!!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_path = './MiniBlogData';
file_format = 'csv';
time_format = 2;
time_scale = 1;
time_offset = 0;

if strcmp(file_format, 'txt')==0 && strcmp(file_format, 'csv')==0
    warning('The format might be unsupported.')
end

list = dir(sprintf('%s/*.%s', file_path, file_format));
if isempty(list)
    warning('None of supported files are found.')
end

Data = [];
num_file = length(list);
tic

if num_file>100
    parfor n = 1:num_file
        tic
        data = readtable(sprintf('%s/%s', file_path, list(n).name));
        Data = [Data; data];
        if mod(n, 100)==0 || n==num_file
            fprintf('File %d/%d, time=%.2fsec\n', n, num_file, toc);
        end
    end

else
    for n = 1:num_file  %num_file = 1;

        data = readtable(sprintf('%s/%s', file_path, list(n).name));
        Data = [Data; data];
        if mod(n, 100)==0 || n==num_file
            fprintf('File %d/%d, time=%.2fsec\n', n, num_file, toc);
        end
    end
end


Seqs = struct('Time', [], ...
              'Mark', [], ...
              'Start', [], ...
              'Stop', [], ...
              'Feature', []);
 
Seqs(1).Start = 0;


Events = Data(:,2);  %Mark
Events = table2cell(Data(:,2));


Time = table2cell(Data(:,1));



for i = 60000:size(Data,1)
    % sequence ID  only 1 sequences is considered here.
    id = 1;
    
    % event ID
    Seqs(id).Mark = Events;
    for k = 1:length(Events);
         Seqs(id).Mark{k} = Events{k};
    end

    %Seqs(id).Mark = str2num(Seqs(id).Mark);
    
    % time stamp
    time = Time{i};
    switch time_format
        case 1 % real number
            t = (time-time_offset)/time_scale;

        case 2 % "year/month/day , XXXX/XX/XX "
            Month = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
            Month = cumsum(Month);
            year = str2double(time(1:4));
            month = str2double(time(6:7));
            day = str2double(time(9:10));
            
            %t = (year-1)*365*24*60 + (Month(month)+day-1)*24*60
            t = (year-1)*365 + (Month(month)+day-1);
            t = (t-time_offset)/time_scale;
        otherwise
            warning('input the right time_format selection.')
    end
    Seqs(id).Time = [Seqs(id).Time; t];

    if mod(i, 100)==0 || i == size(Data,1)
        fprintf('Events %d/%d, time=%.2fsec\n', i, size(Data,1), toc);
    end
end

%L = zeros(length(Seqs),1);
if length(Seqs)>100
    parfor id = 1:length(Seqs)
        Seqs(id).Start = max([Seqs(id).Time(1)-eps, 0]);
        Seqs(id).Stop = Seqs(id).Time(end)+eps;  %end -> last index of a array
    end
else
    for id = 1:length(Seqs)
        Seqs(id).Start = max([Seqs(id).Time(1)-eps, 0]);
        Seqs(id).Stop = Seqs(id).Time(end)+eps;    
    end
end

fprintf('Finish! Time=%.2fsec\n', toc);

%end

