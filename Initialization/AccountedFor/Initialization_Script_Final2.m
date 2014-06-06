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
    clear all
    fprintf('Please run mouseGO.py')
end