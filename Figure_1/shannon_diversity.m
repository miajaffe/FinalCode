%% Shannon Diversity
% This script calculates the Shannon-Weiner indices for each sample as well 
% as the number of unique proteins in each sample. The
% Shannon index takes into account both the richness and the evenness of a
% community. Thus, a sample with a greater number of proteins as well as a
% more equal distribution of protein abudance will be more diverse.

clear all
close all hidden
clc
load('../Initialization/axes.mat');
load('../Initialization/normOverlord_shannon.mat')
normOverlord = normOverlord_shannon * 1000;
%% Reformat matrix so each row is a protein and each column is a sample
all_samples = [];
all_labels = {};
proteins = axes{1}'; %generates list of protein ids
for mouse_num = 1:3
    for colonization = 1:3
        for loc = 1:5
            all_samples = [all_samples normOverlord(:,mouse_num, colonization, loc)];
            label = strcat(axes{2}{mouse_num}, '_', axes{3}{colonization} , '_', axes{4}{loc});           
            all_labels = [all_labels label];
        end
    end
end


%%
% reshape1 matrix looks like this:
% mouse1 coloniz 1 loc 1-5
% mouse1 coloniz 2 loc 1-5
% mouse1 coloniz 3 loc 1-5
% mouse 2 coloniz 1 loc 1-5
% mouse 2 coloniz 2 loc 1-5
% mouse 2 coloniz 3 loc 1-5
% mouse 3 coloniz 1 loc 1-5
% mouse 3 coloniz 2 loc 1-5
% mouse 3 coloniz 3 loc 1-5

[H,VarH]=index_SaW(all_samples); %This helper method taken from MatLab File Exchange calculates the shannon diversity of each sample and returns the indices.
reshape1 = reshape(H, 5,9)'; %9 by 5 where each column is a different region along the GI tract
reshape1 = [reshape1(:,5) reshape1(:,3) reshape1(:,2) reshape1(:,1) reshape1(:,4)] %put in order of GI tract
coloniz_1  = mean(reshape1([1 4 7],:)); %takes average of replicates
sem_coloniz_1 = std(reshape1([1 4 7],:))./sqrt(3); % calculate standard error of the mean

coloniz_2  = mean(reshape1([2 5 8],:));
sem_coloniz_2 = std(reshape1([2 5 8],:))./sqrt(3);

coloniz_3  = mean(reshape1([3 6 9],:));
sem_coloniz_3 = std(reshape1([3 6 9],:))./sqrt(3);
sem_matrix = [sem_coloniz_1' sem_coloniz_2' sem_coloniz_3'];
matrix_for_bar_graph = [coloniz_1' coloniz_2' coloniz_3']; %5 by 3 matrix where 5 is locations and 3 is the colonization states


%% anova - cecum
anova_matrix_cecum = [reshape1([1 4 7],4) reshape1([2 5 8],4) reshape1([3 6 9],4)];
[p_anova, table, stats] = anova1(anova_matrix_cecum);
figure
sig_cecum = multcompare(stats) %Tukey test
%% anova - prox colon
anova_matrix_colon = [reshape1([1 4 7],5) reshape1([2 5 8],5) reshape1([3 6 9],5)];
[p_anova2, table2, stats2] = anova1(anova_matrix_colon);
figure
sig_colon = multcompare(stats2) %Tukey test

%% Plot
colormap('jet')
bar_graph = barweb(matrix_for_bar_graph, sem_matrix);
% set(bar_graph,'FaceColor','r');
legend('Germ Free', 'B. Theta', 'Conventional')
title('Shannon Diversity of Samples')
ylabel('Shannon-Weiner Index')


%% calculate the number of unique proteins in each sample
protein_num= zeros(45,1);
% calculate number of prot in each sample
for sample = 1:45
   non_zero_indices = find(all_samples(:,sample));
   prot_num(sample) = length(non_zero_indices);
end
prot_num_reshaped = reshape(prot_num, 5,9)';
prot_num_reshaped = [prot_num_reshaped(:,5) prot_num_reshaped(:,3) prot_num_reshaped(:,2) prot_num_reshaped(:,1) prot_num_reshaped(:,4)]; %put in order of GI tract
coloniz_1_prot_num  = mean(prot_num_reshaped([1 4 7],:)); %takes average of replicates
sem_coloniz_1_prot = std(prot_num_reshaped([1 4 7],:))./sqrt(3); %calculate stand error of mean
coloniz_2_prot_num  = mean(prot_num_reshaped([2 5 8],:)); %takes average of replicates
sem_coloniz_2_prot = std(prot_num_reshaped([2 5 8],:))./sqrt(3);
coloniz_3_prot_num  = mean(prot_num_reshaped([3 6 9],:)); %takes average of replicates
sem_coloniz_3_prot = std(prot_num_reshaped([3 6 9],:))./sqrt(3);

sem_matrix_prot = [sem_coloniz_1_prot' sem_coloniz_2_prot' sem_coloniz_3_prot'];
protein_num_matrix = [coloniz_1_prot_num' coloniz_2_prot_num' coloniz_3_prot_num'];%5 by 3 matrix where 5 is locations and 3 is the colonization states
figure
bar_graph = barweb(protein_num_matrix, sem_matrix_prot);
legend('Germ Free', 'B. Theta', 'Conventional')
ylabel('Unique Proteins') 

%% anova - cecum - number proteins - GF has  more proteins than conventional
anova_cecum_nprot = [prot_num_reshaped([1 4 7],4) prot_num_reshaped([2 5 8],4) prot_num_reshaped([3 6 9],4)];
[p_anova3, table3, stats3] = anova1(anova_cecum_nprot);
figure
sig_cecum_nprot = multcompare(stats3)
%% anova - prox colon - number proteins - GF has more proteins than conventional
anova_colon_nprot = [prot_num_reshaped([1 4 7],5) prot_num_reshaped([2 5 8],5) prot_num_reshaped([3 6 9],5)];
[p_anova4, table4, stats4] = anova1(anova_colon_nprot);
figure
sig_colon_nprot = multcompare(stats4)

%% anova - stomach - not significant
anova_stomach_nprot = [prot_num_reshaped([1 4 7],1) prot_num_reshaped([2 5 8],1) prot_num_reshaped([3 6 9],1)];
[p_anova5, table5, stats5] = anova1(anova_stomach_nprot);

%% anova - jejunum - not significant
anova_jejunum_nprot = [prot_num_reshaped([1 4 7],2) prot_num_reshaped([2 5 8],2) prot_num_reshaped([3 6 9],2)];
[p_anova6, table6, stats6] = anova1(anova_jejunum_nprot);

%% anova - ileum - not significant
anova_ileum_nprot = [prot_num_reshaped([1 4 7],3) prot_num_reshaped([2 5 8],3) prot_num_reshaped([3 6 9],3)];
[p_anova7, table7, stats7] = anova1(anova_ileum_nprot);

