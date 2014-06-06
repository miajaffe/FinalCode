% This file generates the tab-delimited txt file of the first dimension
% (UniProt ID list) of the axes.mat file for use in mouseGO.py.
load('axes')
fileID = fopen('axes.txt','w');
formatSpec = '%s\t';
for i = 1:1:size(axes{1},2)
    fprintf(fileID,formatSpec,axes{1}{i});
end
fclose(fileID);