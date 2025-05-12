function [task_pre, run_loc] = loadData_v2(subject, path)
dirInfo = dir(path);

task_pre = struct('data', [], 'header', [], 'eof_eeg', [],'behavior', [], 'eof_behavior', []);

dirNames = {dirInfo.name};
dirFolders = {dirInfo.folder};
run_loc = zeros(length(dirNames),1);
for i = 1:length(dirNames) %go through each run
    exName = dirNames{i};
    exFolder = dirFolders{i};
    fileName = [exFolder '/' exName];
    [tempSignal, header] = sload([fileName]);

        if (isempty(task_pre.header))
            task_pre.header = header;
        else
            task_pre.header.EVENT.POS = cat(1, task_pre.header.EVENT.POS, header.EVENT.POS+size(task_pre.data, 1));
            task_pre.header.EVENT.TYP = cat(1, task_pre.header.EVENT.TYP, header.EVENT.TYP);
        end
        task_pre.data = cat(1, task_pre.data, tempSignal);
        task_pre.eof_eeg = cat(1, task_pre.eof_eeg, size(task_pre.data, 1));
        if i == 1 && run_loc(i) == 0
            run_loc(i) = size(tempSignal,1);
        else
            run_loc(i) = run_loc(i-1) + size(tempSignal,1);
        end
end
delete sopen.mat
end
