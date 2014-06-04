% Initialization Script
% Run me!
clear all
close all
clc
filename = 'Longitudinal_RawCounts_Final.xlsx';
[OverlordMatrix,PeptideMap,LetterMap,axes] = PrepareRawData(filename);
save('OverlordMatrix.mat','OverlordMatrix')
save('axes','axes');
save('PeptideMap.mat','PeptideMap');
save('LetterMap.mat','LetterMap');
