function [ output ] = FLP_Fuzzification( FuzzySet, CrispInput )
% FLP_Fuzzification Calculates the fuzzy antecedent membership grades
%
% This function converts each of the crisp input values to a fuzzy
% membership grade for each item in the set. For example, if there are 10
% crisp input values and five items in a set, the resulting output is a 10
% by 5 matrix with the membership grades for each item being in a column
% and the rows being each crisp value. If there are multiple sets in the
% FuzzySet, the output will be a cell array with each cell containing the
% matrix of membership grades for that set.
%
% Input
% FuzzySet - a FuzzySet object generated by FLP_LoadFuzzySets
% CrispInput - a CrispInput object generated by FLP_LoadCrispInput
%
% Output
% output - the antecedent membership grades
%
% Author: Jim Kunce (jdk_acct@yahoo.com)

output = cell(FuzzySet.Count-1,1); % pre-allocate the output

for i = 1:FuzzySet.Count % loop through each set
    
    if ~strcmp(FuzzySet.Set{i,1},'Output') % do not calculate grades for the Output set
    
        mfArray = zeros(size(CrispInput,1),FuzzySet.ItemCount(i,1)); % initialize the membership grade output matrix
        
        for j = 1:FuzzySet.ItemCount(i,1) % loop through each item in the set
            
            mfArray(:,j) = FLP_trapzMF(CrispInput(:,i),FuzzySet.Parms{i,1}(j,:)); % calc the membership grades for each CrispInput
            
        end
    
        output{i,1} = mfArray; % store the grades to a cell array
        
    end
        
end

end

