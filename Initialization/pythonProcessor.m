% Creates a 3 column array, with as many rows as there are proteins as
% listed in element 1 of axes from OverlordMatrix.  Each row in this array
% (GOArray) corresponds to the protein of the same row in axes; note that
% there will be empty rows in GOArray for proteins without GO codes.
clear all
close all
clc
% Opens Python data
filename = fopen('final_Dictionary_GO.txt');
% Reads information from Python file
GoDictionary = textscan(filename,'%s','delimiter','\n');
counter = 0;
counter2 = 0;
% Generates an array corresponding to 'axes' from OverlordMatrix: protein 2
% has the GO codes in row 2 of GOArray and so on. Column 1 has the string
% form of the GO ID's, column 2 are the associated labels, and column 3 has
% the integer form of the GO IDs
for ii = 1:1:size(GoDictionary{1},1)
    counter = counter + 1;
    % Sets the current row in GOArray, corresponding to a protein in 'axes'
    if counter == 1
        currentIndex = str2num(GoDictionary{1}{ii});
    elseif counter == 2
        tempValue = GoDictionary{1}{ii};
        tempValue = tempValue(2:length(tempValue) - 1);
        for iii = 1:1:((length(tempValue) + 2) / 14)
            % Generates a cell array of the GO code ID's, e.g. 'GO:0000001'
            value{(ii-2)/4 + 1,1}{iii} = tempValue(((iii-1)*14) + 2: ((iii-1)*14) + 11);
            % Generates a cell array of the numerical GO code ID's, e.g.
            % 'GO:0000001' to 1
            value{(ii-2)/4 + 1,3}(iii) = str2num(tempValue(((iii-1)*14) + 5: ((iii-1)*14) + 11));
        end
    elseif counter == 3
        tempValue = GoDictionary{1}{ii};
        tempValue = tempValue(2:length(tempValue) - 1);
        for iii = 1:1:((length(tempValue) + 2) / 5)
            % Generates a cell array of the GO code labels (F=molecular function, 
            % C=cellular component, P=biological process)
            value{(ii-3)/4 + 1,2}{iii} = tempValue(((iii-1)*5) + 2);
        end
    else
        % Fills in the current row of the GOArray with the string GO ID,
        % the labels, and the numerical GO ID
        counter2 = counter2 + 1;
        counter = 0;
        GOArray{currentIndex,1} = value{counter2,1};
        GOArray{currentIndex,2} = value{counter2,2};
        GOArray{currentIndex,3} = value{counter2,3};
    end
end
% Saves resuting GOArray to file
save('GOArray.mat')