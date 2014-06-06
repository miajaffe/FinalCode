%% Summary Stats - Pie Charts
% This file generates the pie charts that describe the percent abundance of
% each protein. The percent abundance of each protein is calculated, and
% all proteins with < 2% abundance are grouped into an "other" category in
% the pie chart. 

clear all
close all hidden
clc
load('axes.mat');
load('normOverlordFinal.mat'); 
location_labels = axes{4};
colonization_labels = axes{3};
normOverlord = normOverlordFinal * 100; 
proteins = axes{1}' %list of all proteins

final_table = {};
final_names = {};
count = 1;      
cutoff = 0.02; %cutoff protein abundance to group into "other" category
% Each protein was manually tied to a specific color. Following are the
% list of colors used in the pie charts (N.B. Some colors may have been
% changed using a GUI for our final figures).
colors =[166 206 227; 31 120 180; 178 223 138; 51 160 44; 251,154,153; 227,26,28; 253,191,111; 255,127,0; 202,178,214; 106,61,154; 255,255,153; 177,89,40];
colors = colors./255;
colors2 = [27,158,119; 217,95,2; 117,112,179; 231,41,138; 102,166,30; 230,171,2; 166,118,29; 0 0 0.75; 0 0.75 0; 0.75 0 0; 1 1 0; 1 0 1; 0 1 1; 0 0 0.5; 0 0.5 0; 0.5 0 0; 0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0];
colors2 = colors2./255;

abundances ={};
for location = 1:5
    for colonization = 1:3
        title_label = strcat(location_labels(location), colonization_labels(colonization));
        % sum up protein abundance for all of the replicates
        replicate_sum = sum(normOverlord(:,:,colonization,location), 2); 
        % total protein abundance for all proteins 
        total_prot_abundance = sum(replicate_sum, 1);
        n_proteins = length(replicate_sum);
        total_abundance_vect = repmat(total_prot_abundance, n_proteins, 1);
        %calculate the percentage of each protein abundance compared to the
        %total
        prot_percentages = replicate_sum ./ total_abundance_vect;
        % Group all of the proteins with < 2% abundance into one other
        % group and add up total percentage of the group
        [rows, cols] = find(prot_percentages <= cutoff);
        prot_small_abundance = replicate_sum(rows,1);
        summed_small_abundance = sum(prot_small_abundance);
        % Group all of the proteins with >2% abudances into a matrix
        [rows2, cols2] = find(prot_percentages > cutoff);
        fprintf('%d', count);
        %Generate a list of the associated protein names
        large_prot = proteins(rows2,1)
        prot_large_abundance = replicate_sum(rows2,1)
        prot_name_counts = {large_prot prot_large_abundance} %names and corresponding counts
        final_table = [final_table prot_name_counts]; %The final table includes the names/counts for each sample
        final_abundances = [prot_large_abundance; summed_small_abundance];
        abundances = [abundances final_abundances];
        %figure 
        subplot(5,3,count)
        abundance = pie(final_abundances) % create pie charts
        title(title_label(1));
        count = count + 1;
        
        %set colors of the pie chart depending on the protein in the chart        
        hp = findobj(abundance, 'Type', 'patch');
        for i = 1: length(large_prot)
            switch(large_prot{i})
                case 'A8DUK4'
                    set(hp(i), 'FaceColor', colors(1,:));
                case 'B2RS76'
                    set(hp(i), 'FaceColor', colors(2,:));
                case 'P02815'
                    set(hp(i), 'FaceColor', colors(3,:));
                case 'P05208'
                    set(hp(i), 'FaceColor', colors(4,:));
                case 'P07146'
                    set(hp(i), 'FaceColor', colors(5,:));
                case 'P07743'
                    set(hp(i), 'FaceColor', colors(6,:));
                case 'Q3SYP2'
                    set(hp(i), 'FaceColor', colors(7,:));
                case 'Q61900'
                    set(hp(i), 'FaceColor', colors(8,:));
                case 'Q62395'
                    set(hp(i), 'FaceColor', colors(9,:));
                case 'Q64097'
                    set(hp(i), 'FaceColor', colors(10,:));
                case 'Q6P8U6'
                    set(hp(i), 'FaceColor', colors(11,:));
                case 'Q7TNM8'
                    set(hp(i), 'FaceColor', colors(12,:));
                case 'Q7TPZ8'
                    set(hp(i), 'FaceColor', colors2(1,:));
                case 'Q8C5B4'
                    set(hp(i), 'FaceColor', colors2(2,:));
                case 'Q8K0C5'
                    set(hp(i), 'FaceColor', colors2(3,:));
                case 'Q8R0I0'
                    set(hp(i), 'FaceColor', colors2(4,:));
                case 'Q91WB5'
                    set(hp(i), 'FaceColor', colors2(5,:));
                case 'Q91WL7'
                    set(hp(i), 'FaceColor', colors2(6,:));
                case 'Q91X79'
                    set(hp(i), 'FaceColor', colors2(7,:));
                case'Q91XA9'
                    set(hp(i), 'FaceColor', colors2(8,:));
                case 'Q921Y7'
                    set(hp(i), 'FaceColor', colors2(9,:));
                case 'Q9CPN9'
                    set(hp(i), 'FaceColor', colors2(10,:)); 
                case 'Q9CQ52'
                    set(hp(i), 'FaceColor', colors2(11,:));
                case 'Q9CQC2'
                    set(hp(i), 'FaceColor', colors2(12,:));
                case 'Q9CR35'
                    set(hp(i), 'FaceColor', colors2(13,:));
                case 'Q9CY06'
                    set(hp(i), 'FaceColor', colors2(14,:));
                case 'Q9ER05'
                    set(hp(i), 'FaceColor', colors2(15,:));
                case 'Q9JM84'
                    set(hp(i), 'FaceColor', colors2(16,:));
                case 'Q9R0T7'
                    set(hp(i), 'FaceColor', colors2(17,:));
                case 'Q9CQC2'  
                    set(hp(i), 'FaceColor', colors2(18,:));
                case 'P02816'
                    set(hp(i), 'FaceColor', colors2(19,:));    
            end
        end
    end
end
      