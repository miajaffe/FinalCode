% Initialization Script
% Run me! This script generates the first half of helper files required for
% all downstream scripts in the Figure_X folders, where X is an integer.
% Please see the documentation for the individual scripts for a detailed
% description of each step. It is required to run mouseGO.py right after
% this script.
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
[normOverlord_shannon,redundancies_shannon] = overlord_normalizer_supply_matrix();
clear all
run('axesTextGenerator.m');
clear all
clc
% It is now necessary to run mouseGO.py