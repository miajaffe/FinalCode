% GET ALL THE GO CODE TITLES. Parses the European Bioinformatics Institute
% site for the GO title associated with each and every occuring GO ID and
% records each GO ID: GO title pair in a Map container named allGODic.
clear all
clc
load GOtoIndexConverter
allnumGO = keys(GOtoIndexConverter);
% Convert numerical GO values to searchable strings
allstrGO = {};
for i = 1:1:length(allnumGO)
    currStr = num2str(allnumGO{i});
    front = 'GO:';
    zeroslen = 7 - length(currStr);
    for j = 1:1:zeroslen
        front = strcat(front,'0');
    end
    allstrGO{i} = strcat(front,currStr);
end

% Store annotation and definition in a Map with GO ID as key
allGODic = containers.Map();
% For each and every GO ID, conduct a url search to retrieve the
% corresponding GO title and add to allGODic with GO ID (string) as the key 
% and the corresponding GO title as the value.
for i = 1:1:size(allstrGO,2)
    key = allstrGO{i};
    url = strcat('http://www.ebi.ac.uk/QuickGO/GTerm?id=',key);
    tempData = urlread(url);
    startI = strfind(tempData,'<td class="label">Definition') + 40;
    counter = startI;
    definition = '';
    while tempData(counter) ~= '<';
        definition = strcat(definition,tempData(counter));
        counter = counter + 1;
    end
    startI = strfind(tempData,'<title>') + 18;
    annot = '';
    counter = startI;
    while tempData(counter) ~= '<';
        annot = strcat(annot,tempData(counter));
        counter = counter + 1;
    end
    allGODic(key) = {annot,definition};
    fprintf('%d\n',i)
end
save('allGODic.mat','allGODic')