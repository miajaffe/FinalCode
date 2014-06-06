
%% Figure 2
% Hierarchical Clustering

% Figure 2.B.
% This script runs a Hierarchical Clustering Analysis on the proteomic 
% database and generations a clustergram.

% Load Initialized Data
clear all;
close all hidden;
clc;
load('../Initialization/axes.mat');
load('../Initialization/normOverlordFinal.mat');
all_samples = [];
all_labels = {};

% The following code creates a 2D matrix of proteins vs samples and creates
% a clustergram to visualize hierarchical clustering of the data.

% This creates a clustergram of the data where proteins are on the y axis
% and the samples are on the x axis. The goal is to group samples into
% individual clusters. 


% Loops through all mouse replicates, colonization states, and locations to
% create a 2D matrix where each row is a protein_id and each column is a
% different sample (sample A, B, etc.). The matrix, all_samples, contains 
% the abundance of each protein at each sample.  
% The all_labels matrix is a 1x45 matrix
% that contains labels for each sample in the order that they are in the
% all_samples matrix. The all_labels matrix and the proteins matrix can be 
% used for labeling the axes of any plots.

for mouse_num = 1:3;
    for colonization = 1:3;
        for loc = 1:5;
            all_samples = [all_samples normOverlordFinal(:,mouse_num, colonization, loc)];
            label = strcat(axes{2}{mouse_num}, '_', axes{3}{colonization} , '_', axes{4}{loc});           
            all_labels = [all_labels label];
        end;
    end;
end;
      
proteins = axes{1}'; %generates list of protein ids

plot(clustergram(all_samples, 'RowLabels', proteins, 'ColumnLabels', all_labels', 'DisplayRange', 0.0015, 'Symmetric', 'true', 'colormap', 'winter'));
colorbar

print -dpdf -r600 Figure_2B.pdf


