function [ consequentGrades ] = FLP_FuzzyRuleEval( AntMemberGrades, FuzzyRules )
% FLP_FuzzyRuleEval Evaluates Fuzzy Rules & Consequent Membership Grades
%
% Based on the antecedent membership grades, each of the fuzzy rules is
% evaluated using the Mamdani procedure. The Mamdani procedure consists of
% taking the minimum antecedant grade as the membership grade for the rule.
% The consequent membership grade is then determined by taking the maximum
% of the rule grades for each consequent.
%
% Input
% AntMemberGrades - a antecedent membership grade object generated by FLP_Fuzzification
% FuzzyRules - a fuzzy rules object generated by FLP_LoadFuzzyRules
%
% Output
% the consequent membership grades
%
% Author: Jim Kunce (jdk_acct@yahoo.com)

setCt = size(AntMemberGrades,1); % get the number of fuzzy sets
inputCt = size(AntMemberGrades{1,1},1); % get the number of crisp inputs
ruleCt = size(FuzzyRules,1); % get the number of rules
maxLevelOutput = max(FuzzyRules(:,end)); % get the number of enumerated output levels
consequentGrades = zeros(inputCt,maxLevelOutput); % pre-allocate the output

for i = 1:inputCt % loop thru each crisp input
    
    membershipGrades = zeros(ruleCt,setCt); % pre-allocate a matrix to hold the antecedent grades
    
    for j = 1:setCt % loop thru each set
       
        levelCt = size(AntMemberGrades{j,1}(i,:),2); % get the number of levels in the set
        levelGrades = repmat(AntMemberGrades{j,1}(i,:),ruleCt,1); % replicate the level grades row for each rule
        
        levelMap = repmat(1:levelCt,ruleCt,1); % create a row of enumerated level #'s for each rule 
        levelLookup = repmat(FuzzyRules(:,j),1,levelCt); % replicate level # for rule for each level
        levelLI = levelLookup==levelMap; % create a logical index matching the level # for the rule
        membershipGrades(:,j) = max(levelGrades.*levelLI,[],2); % extract the level grades for the rule
        membershipGrades(FuzzyRules(:,j)==0,j) = NaN; % insert NaN when no level for the set is present in the rule
        
    end
    
    ruleGrades = repmat(min(membershipGrades,[],2),1,maxLevelOutput); % take the minimum of the antecedents
    outputClassMap = repmat(1:maxLevelOutput,ruleCt,1); % create a role of enumerated output levels for each rule
    outputClassLookup = repmat(FuzzyRules(:,end),1,maxLevelOutput); % replicate the output level for the rule for each output level
    outputClassLI = outputClassMap==outputClassLookup; % create a logical index matching the output level # for the rule
    consequentGrades(i,:) = max(ruleGrades.*outputClassLI,[],1); % take the max of the grades as the consequent membership grade
    
end

end

