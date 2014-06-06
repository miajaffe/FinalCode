% Initialization Script 2
% Run me! This script generates the second half of helper files required for
% all downstream scripts in the Figure_X folders, where X is an integer.
% Please see the documentation for the individual scripts for a detailed
% description of each step. The output file of mouseGO.py
% (final_Dictionary_GO.txt) is required to be in the same directory for
% this code to run.
if exist('final_Dictionary_GO.txt', 'file') > 0
    clear all
    run('pythonProcessor.m')
    clear all
    run('GOTally.m')
    clear all
    run('GOTally_shannon.m')
    clear all
    run('getAllGODef.m')
    clear all
else
    % Error message if final_Dictionary_GO.txt not found
    clear all
    fprintf('Please run mouseGO.py')
end