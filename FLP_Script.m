% Script to run Fuzzy Logic Processor example
%
% Author: Jim Kunce (jdk_acct@yahoo.com)

FuzzySet = FLP_LoadFuzzySets('WashingMachine_FuzzySets.csv');
CrispInput = FLP_LoadCrispInput(FuzzySet,'WashingMachine_CrispInputs.csv');
FuzzyRules = FLP_LoadFuzzyRules(FuzzySet,'WashingMachine_FuzzyRules.csv');

AntMemberGrades = FLP_Fuzzification(FuzzySet, CrispInput);
ConsqMemberGrades = FLP_FuzzyRuleEval(AntMemberGrades,FuzzyRules);
[CrispOutput, OutputMF, X] = FLP_DeFuzzification(ConsqMemberGrades, FuzzySet, 100);
plotdata = FLP_3DMFplot( FuzzySet, FuzzyRules, 1, 2, 20 );