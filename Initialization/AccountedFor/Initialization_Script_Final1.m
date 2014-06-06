% Initialization Script
% Run me!
clear all
close all
clc
filename = 'Longitudinal_RawCounts_Final.xlsx';
[OverlordMatrix,PeptideMap,LetterMap,axes] = PrepareRawData(filename);
save('OverlordMatrix.mat','OverlordMatrix')
save('axes.mat','axes');
save('PeptideMap.mat','PeptideMap');
save('LetterMap.mat','LetterMap');
clear all
[normOverlord,redundancies] = OverlordNormalizer('MusProtRaw.fasta');
clear all
[normOverlord_shannon,redundancies] = overlord_normalizer_supply_matrix();
clear all
run('axesTextGenerator.m');
clear all
clc
% It is now necessary to run mouseGO.py