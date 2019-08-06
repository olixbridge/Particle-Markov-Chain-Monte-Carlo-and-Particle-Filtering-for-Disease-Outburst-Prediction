// compile with R CMD SHLIB  testMMHPMCMCSkeletonForDifferentialEquations.c

#include <assert.h>

void integrateODEFromPreviousToCurrentTimepointModelSpecificFullyQualified(int iTimepoint, 
                     int countParameters, 
                     double parameters[countParameters],
                     int countStateVariables, 
                     int countObservationTimepoints,
                      double arrayTimesForTimepoints[countObservationTimepoints],
                      double integrationTimestepParm,
                      // the below is both an INPUT and an OUTPUT variable!
            double stateVariablesForTimepoint[countObservationTimepoints][countStateVariables]);

extern void initializeStateVariablesAtIndexToUniformValue(int iStartingIndexToInitialize,
                                     int countValuesToInitialize,
                                     double initialValueToImpose,
                                     int countStateVariables,
                                     double stateVariablesForParticleSampled[countStateVariables]);

extern void initializeStateVariablesAtIndexToSubstateSpecificValues(int iStartingIndexToInitialize,
                                     int countValuesToInitialize,
                                     double arrayInitialValuesToImpose[countValuesToInitialize],
                                     int countStateVariables,
                                     double stateVariablesForParticleSampled[countStateVariables]);

#define initial_remandPopulationFemale_child_min (30)
#define initial_remandPopulationFemale_teen_min (30)
#define initial_remandPopulationFemale_old_min (30)
#define initial_remandPopulationFemale_child_max (40)
#define initial_remandPopulationFemale_teen_max (40)
#define initial_remandPopulationFemale_old_max (40)


#define initial_SentencedPopulationFemale_child_min (50)
#define initial_SentencedPopulationFemale_teen_min (50)
#define initial_SentencedPopulationFemale_old_min (50)
#define initial_SentencedPopulationFemale_child_max (100)
#define initial_SentencedPopulationFemale_teen_max (100)
#define initial_SentencedPopulationFemale_old_max (100)

#define initial_remandPopulationMale_child_min (500)
#define initial_remandPopulationMale_teen_min (500)
#define initial_remandPopulationMale_old_min (500)
#define initial_remandPopulationMale_child_max (600)
#define initial_remandPopulationMale_teen_max (600)
#define initial_remandPopulationMale_old_max (600)

#define initial_SentencedPopulationMale_child_min (700)
#define initial_SentencedPopulationMale_teen_min (700)
#define initial_SentencedPopulationMale_old_min (700)
#define initial_SentencedPopulationMale_child_max (1000)
#define initial_SentencedPopulationMale_teen_max (1000)
#define initial_SentencedPopulationMale_old_max (1000)


#define Month (30.0)
#define Year (365.0)
#define Week (7.0)

#define TotalPopulation (1104825)
#define FractionOfMale (0.49)

// Female (0), Male(1)
#define Dimensions_Sex (2)
#define Dimensions_Age (3)
#define Female_child_index (0)
#define Female_teen_index (1)
#define Female_old_index (2)
#define Male_child_index (0)
#define Male_teen_index (1)
#define Male_old_index (2)

// define the start index of each stock
#define GeneralPopulation_StateIndex (0)
#define InitialDetentionPopulation_StateIndex (6)
#define RemandPopulation_StateIndex (12)
#define SentencedCustody_StateIndex (18)
#define OnBail_StateIndex (24)
#define OnProbation_StateIndex (30)

// define the index of each stocks
//need to change to 36 instead of 12
#define GeneralPopulation_Female_StateIndex (0)
#define GeneralPopulation_Male_StateIndex (1)
#define InitialDetentionPopulation_Female_StateIndex (2)
#define InitialDetentionPopulation_Male_StateIndex (3)
#define RemandPopulation_Female_StateIndex (4)
#define RemandPopulation_Male_StateIndex (5)
#define SentencedCustody_Female_StateIndex (6)
#define SentencedCustody_Male_StateIndex (7)
#define OnBail_Female_StateIndex (8)
#define OnBail_Male_StateIndex (9)
#define OnProbation_Female_StateIndex (10)
#define OnProbation_Male_StateIndex (11)

#define integrationTimestep (0.01)

// dynamic parameters
//dont change
#define log_arrestHazard_PF_SEIRStateIndex  (12) 
#define logit_fractionOfInitallyDetainedEnteringRemand_PF_StateIndex (13)

#define arrestHazardRate_PFMin (0.005) // male
#define arrestHazardRate_PFMax (0.01) // male 
#define arrestHazardRate_fractionOfFemaleToMale (0.1)
#define arrestHazard_RandomWalkStdDev (0.8) //0.4

#define fractionOfInitallyDetainedEnteringRemand_PFMin (0)
#define fractionOfInitallyDetainedEnteringRemand_PFMax (0.5)
#define fractionOfInitallyDetainedEnteringRemand_RandomWalkStdDev (0.3) //0.4

// MCMC sampled parameters
#define fractionOfInitiallyDetainedReturningImmediatelyToCommunityParaMin (0.1)
#define fractionOfInitiallyDetainedReturningImmediatelyToCommunityParaMax (0.9)

#define logit_fractionOfInitiallyDetainedReturningImmediatelyToCommunity_parameterIndex (0)


//have 3 extra, have three
#define Empirical_Remand_Male_Index (0)
#define Empirical_Remand_Female_Index (1)
#define Empirical_SentencedCustody_Male_Index (2)
#define Empirical_SentencedCustody_Female_Index (3)
