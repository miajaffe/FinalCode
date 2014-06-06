% Examines to find GO titles for which there is a statistically significant
% difference in relative abundance between any pair within the three 
% colonization states. This is done by first running anova1 and
% subsequently running mafdr and selecting GO titles below a q-treshold of
% 0.01.
clear all
close all
% clc
load('../Initialization/GOenrichMat')
load('../Initialization/GOenrichMat_shannon')
load('../Initialization/axes')
load('../Initialization/GOtoIndexConverterStr')
load('../Initialization/IndextoGOConverterStr')
load('../Initialization/allGODic')
load('../Initialization/titleToGO')
%% Find GO codes with q-value <= 0.01
allGO = keys(GOtoIndexConverterStr);
% Computes the p-value associated with anova1 between all three
% colonization states.  Each colonization state is considered to be a
% mutually independent observation (column) and each colonization state has
% 3 dependent observations (rows/mouse replicates). Note that this loop
% first splits GOenrichMat_shannon into its 3 colonization states then sums
% across all locations.
for i = 1:1:length(allGO)
    GOcurr = allGO{i};
    for j = 1:1:3
        index = GOtoIndexConverterStr(GOcurr);
        GF_BT_RF_GO1(j,1) = sum(GOenrichMat_shannon(index,j,1,:),4);
        GF_BT_RF_GO1(j,2) = sum(GOenrichMat_shannon(index,j,2,:),4);
        GF_BT_RF_GO1(j,3) = sum(GOenrichMat_shannon(index,j,3,:),4);
    end
    [p,table,stats] = anova1(GF_BT_RF_GO1,{'GF','BT','RF'},'off');
    pall(i) = p;
    statsall(i) = stats;
end
% Sort vector of all anova1 p-values
pall_sorted = sort(pall,'ascend');
% Generate associated false discovery rate and q-value of the anova1
% p-values
[FDR, q] = mafdr(pall);
[sortedq, qind] = sort(q,'ascend');
counter = 1;
% q-value threshold
qcutoff = 0.01;
statssorted = statsall(qind);
% Checks each GO title (sorted by the respective q-values in an ascending
% fashion) to see if its q-values are below the q-value cutoff.  If so,
% records the GO title, the means and stds for the abundances of each
% colonization state, the anova1 p-value, and the q-value.
while sortedq(counter) <= qcutoff
    tempnum = mod(qind(counter),3);
    index = qind(counter);
    value = allGO(index);
    temp = allGODic(value{1});
    counts = titleToGO(temp{1});
    counts = counts{3};
    significantGO{counter,1} = temp{1};
    significantGO{counter,2} = mean(sum(counts(1,:,1,:),4));
    significantGO{counter,3} = mean(sum(counts(1,:,2,:),4));
    significantGO{counter,4} = mean(sum(counts(1,:,3,:),4));
    significantGO{counter,5} = std(sum(counts(1,:,1,:),4));
    significantGO{counter,6} = std(sum(counts(1,:,2,:),4));
    significantGO{counter,7} = std(sum(counts(1,:,3,:),4));
    significantGO{counter,8} = pall_sorted(counter);
    significantGO{counter,9} = sortedq(counter);
    allmeans(counter) = mean([mean(sum(counts(1,:,1,:),4)),mean(sum(counts(1,:,2,:),4)),mean(sum(counts(1,:,3,:),4))]);
    counter = counter + 1;
end
% Sorts GO titles below the threshold according to the mean of their mean
% abundances in a descending fashion.
[sortedallmeans, meanind] = sort(allmeans,'descend');
statsfinal = statssorted(meanind);
% Writes data to file
fileID = fopen('significantGOtitlesQvalsCorrected1pctFDR.csv','w');
formatSpec0 = '%s , %s , %s , %s, %s, %s, %s, %s, %s\n';
header = {'GO Title','GF count','BT count','RF count','GF std','BT std','RF std','p-val','q-val'};
fprintf(fileID,formatSpec0,header{1,:});
formatSpec = '%s,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f\n';
[nrows, ncols] = size(significantGO);
for row = 1:nrows
    fprintf(fileID,formatSpec,significantGO{row,:});
end
fclose(fileID);
% Performs post-processing of each case with q <= q_cutoff to determine
% which pair(s) of colonization states have significantly different
% abundances for the top 15 GO titles (measured by mean abundance).
% Requires use of the GUI to see if each GO code and each colonization
% state thereof is significantly different.
for j = 1:1:15
    figure
    multcompare(statsfinal(j),'alpha',0.01);
end
