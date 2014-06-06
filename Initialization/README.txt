LIST OF REQUIRED FILES:
Initialization_Script_Final1.m
PrepareRawData.m
OverlordMatrixIndexGenerator.m
Longitudinal_RawCounts_Final.xlsx
OverlordNormalizer.m
overlord_normalizer_supply_matrix.m
Initialization_Script_Final2.m
axesTextGenerator.m
mouseGO.py
pythonProcessor.m
GOTally.m
GOTally_shannon.m
getAllGODef.m
GOtitleToGOcode.m
---
System Requirements:
Matlab 2014a
python 2.7.5
---
How to use:
1) Run Initialization_Script_Final1.m in Matlab
2) Run mouseGO.py
3) Run Initialization_Script_Final2.m in Matlab
--
Description of scripts:
Initialization_Script_Final1.m sequentially:
SCRIPT NAME/INPUT(S)/OUTPUT(S) (saved to file)
1) PrepareRawData.m/Longitudinal_RawCounts_Final.xlsx/OverlordMatrix.mat
axes.mat, PeptideMap.mat, LetterMap.mat
2) OverlordNormalizer.m/MusProtRaw.fasta, OverlordMatrix.mat, axes.mat/
normOverlordFinal.mat,redundancies.mat
3) overlord_normalizer_supply_matrix.m/MusProt.mat, OverlordMatrix.mat, axes.mat/
normOverlord_shannon.mat
4) axesTextGenerator.m/axes.mat/axes.txt

It is then required to run mouseGO.py

1) mouseGO.py/gene_association.goa_mouse, axes.txt/final_Dictionary_GO.txt

It is then required to run Initialization_Script_Final2.m
1) pythonProcessor.m/final_Dictionary_GO.txt/GOArray.mat
2) GOTally.m/GOArray.mat, normOverlordFinal.mat, axes.mat/GOenrichMat.mat,
GOtoIndexConverter.mat, GOtoIndexConverterStr.mat, IndextoGOConverter.mat,
IndextoGOConverterStr.mat
3) GOTally_shannon.m/GOArray.mat, normOverlord_shannon.mat, axes.mat/GOenrichMat_shannon.mat
4) getAllGODef.m/GOtoIndexConverter/allGODic.mat
5) GOtitleToGOcode.m/GOenrichMat.mat, axes.mat, GOtindexConverterStr.mat, allGODic.mat/titleToGO.mat
