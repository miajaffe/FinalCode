% Figure S3
% Hierarchical Clustering of GO codes

% Figure S3
% This script runs a Hierarchical Clustering Analysis on the GO Code 
% database and generations a clustergram.


clear all
close all hidden
clc
load('../Initialization/axes.mat');
load('../Initialization/normOverlordFinal.mat');
load('../Initialization/GOenrichMat.mat'); % 2991 * 3 * 3 * 5 4d matrix 
all_samples = [];
all_labels = {};
normOverlord = normOverlordFinal; 

% This loops through all mouse #s, colonization states, and locations to
% create a 2D matrix where each row is a gocode_id and each column is a
% different sample (sample A, B, etc.). The matrix, all_samples, contains the abundance
% of each GO code at each sample.  The all_labels matrix is a 1x45 matrix
% that contains labels for each sample in the order that they are in teh
% all_samples matrix. 
for mouse_num = 1:3
    for colonization = 1:3
        for loc = 1:5
            all_samples = [all_samples GOenrichMat(:,mouse_num, colonization, loc)];
            label = strcat(axes{2}{mouse_num}, '_', axes{3}{colonization} , '_', axes{4}{loc});           
            all_labels = [all_labels label];
        end
    end
end

%simple stats - needed to determine cutoff of clustergram
num_samples = length(GOenrichMat(:,1,1,1)) * 3 * 3 * 5;
reshaped_samples = reshape(all_samples, num_samples,1);
percentages = prctile(reshaped_samples, [25 50 75 90 95])
% percentages =

%         0         0    0.0004    0.0025    0.0063
percentages2 = prctile(reshaped_samples, [81 82 83 84 85])
% percentages2 =

%    0.0008    0.0009    0.0010    0.0011    0.0013

%%
%This clustergram is capturing between 83-84% of the data. In other words
%16-17% of the data is above the 0.0015 cutoff.
%%0.0015
plot(clustergram(all_samples, 'ColumnLabels', all_labels', 'DisplayRange', 0.0025, 'Symmetric', 'true', 'Colormap', winter));

print -dpdf -r600 Figure_S3.pdf
