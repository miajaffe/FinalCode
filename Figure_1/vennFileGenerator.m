%%Prints out tab delimited txt file of the proteins present in the different 
%locations in the GI tract

% Initialization
load axes.mat, load normOverlord_shannon

protIds=axes{1};

% initialize 
% 1 = cecum 2=ileum 3=jejunum 4=proxCol 5=stomach

%find unique protIds from all replicates in each location
indicesOF1ColonizationState=cell(5,3);

uniqueInd=cell(1,5); % unique Indices of proteins of all replicates in a particular colonization 
%state and GI location. For instance the unique indices of all replicates
%in the cecum in Germ free colonization state
uniqueProtIds=cell(3,5); % unique protIds of all replicates in every colonization state
% and GI location
total_indices = {}; %keep tract of all indices to remove across all locations
for GIlocation=1:5
    for colonisationState=1:3
        for mouseReplicate=1:3
            protAbundanceInReplicate=normOverlord_shannon(:,mouseReplicate,colonisationState,GIlocation);
            indicesOF1ColonizationState{GIlocation,mouseReplicate}=find(protAbundanceInReplicate);   
        end
        %saves unique indices of all mouse replicates in a given
        %colonization state in a cell that will contain all unique
        concatInd=cat(1,indicesOF1ColonizationState{GIlocation,1},indicesOF1ColonizationState{GIlocation,2},indicesOF1ColonizationState{GIlocation,3});
        uniqueInd{1,    GIlocation}=unique(concatInd);
        
        %finds associated protIds of the unique indices in a given location
        Id=cell(1,length(uniqueInd{1,GIlocation}));
        indices=uniqueInd{1,GIlocation}; % unique indices of protein Ids
        % Matching unique indices with prot Ids
        for i=1:length(indices)
   
            Id{1,i}=protIds{1,indices(i,1)};
            
        end
        uniqueProtIds{colonisationState,GIlocation}=Id;
        
    end
    total_indices = [total_indices uniqueInd{1,:}]; 
    
    
end
% The total intersection represents the number of proteins found across all
% samples in at least one of the replicates
total_intersection = mintersect(total_indices{1}, total_indices{2}, total_indices{3}, total_indices{4}, total_indices{5}, total_indices{6}, total_indices{7}, total_indices{8}, total_indices{9}, total_indices{10}, total_indices{11}, total_indices{12}, total_indices{13}, total_indices{14}, total_indices{15})
proteins = protIds';
intersect_prot = proteins(total_intersection);
