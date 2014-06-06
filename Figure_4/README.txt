LIST OF REQUIRED FILES:
Figure4.m
getGOcodes.m
GOtitleToGOcode.m
---
System Requirements:
Matlab 2014a
---
How to use:
1) Run Figure4.m in Matlab

--
Description of scripts:
Figure4.m sequentially:
SCRIPT NAME/INPUT(S)/OUTPUT(S) (saved to file)
1) Figure4.m/GOenrichMat.mat, GOenrichMat_shannon.mat, axes.mat, GOtoIndexConverterStr.mat,
IndextoGOConverterStr.mat, allGODic.mat, titleToGO.mat/significantGOtitlesQvalsCorrected1pctFDR.csv

WARNING: There are occasionally GO titles which include commas in their names. It is necessary to 
manually check the resulting csv files for any rows which have an irregular (too many) number of columns.
To fix, merge the GO title names in the excess columns and shift columns to the left.