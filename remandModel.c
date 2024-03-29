// compile with R CMD SHLIB  OsgoodLiuPF2014Model.c
//  This is done via the makefile included.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#undef NDEBUG           // enables assertions

#include <time.h>
//#define MATHLIB_STANDALONE
#include <S.h>
#include <Rdefines.h>
#ifndef USING_R
//#extern void F77_NAME(dqrdca)();
#else
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/Boolean.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <R_ext/RStartup.h>
#endif /* USING_R */

#define DEBUG_LIKELIHOOD
//#define DEBUG_flow
//#define DEBUG_LIKE
//#define DEBUG
#define DEBUG_D
#define DEBUG_stock
// PARAM is only used to test whether it is identity with the anylogic model
//#define PARAM

#include "logging.h"

//extern double min(double a, double b);
double	min(double a, double b);

double	min(double a, double b)
{
	if (a <= b)
		return a;
	else
		return b;
}

#include "remandModel.h"

void initializeStateVariablesAtIndexToUniformValue(int iStartingIndexToInitialize,
                                     int countValuesToInitialize,
                                     double initialValueToImpose,
                                     int countStateVariables,
                                     double stateVariablesForParticleSampled[countStateVariables])
{
  for (int iValueToInitialize = iStartingIndexToInitialize; iValueToInitialize < (iStartingIndexToInitialize + countValuesToInitialize); iValueToInitialize++)
    stateVariablesForParticleSampled[iValueToInitialize] = initialValueToImpose;
}

void initializeStateVariablesAtIndexToSubstateSpecificValues(int iStartingIndexToInitialize,
                                     int countValuesToInitialize,
                                     double arrayInitialValuesToImpose[countValuesToInitialize],
                                     int countStateVariables,
                                     double stateVariablesForParticleSampled[countStateVariables])
{
  for (int iSubstateToInitialize = 0; iSubstateToInitialize < countValuesToInitialize; iSubstateToInitialize++)
    stateVariablesForParticleSampled[iStartingIndexToInitialize+iSubstateToInitialize] = arrayInitialValuesToImpose[iSubstateToInitialize];
}


double resampleInitialParametersForZeroPosteriorModelSpecific(int countParameters,
                    double parameters[countParameters])
{
    assert(countParameters == 1);
    double runif_result_0 = runif(fractionOfInitiallyDetainedReturningImmediatelyToCommunityParaMin, fractionOfInitiallyDetainedReturningImmediatelyToCommunityParaMax);
    parameters[logit_fractionOfInitiallyDetainedReturningImmediatelyToCommunity_parameterIndex] = log(runif_result_0/(1-runif_result_0));
}

double sampleStateVariableInitialStateModelSpecific(int countStateVariables,
                                                    double stateVariablesForParticleSampled[countStateVariables])
{
    assert(countStateVariables == 14);     
    double arrayInitialValueForGeneralPopulation[Dimensions_Sex * Dimensions_Age] = {TotalPopulation * (1.0-FractionOfMale)*0.12, TotalPopulation * FractionOfMale*0.12, TotalPopulation * (1.0-FractionOfMale)*0.06,TotalPopulation * FractionOfMale*0.06, TotalPopulation * (1.0-FractionOfMale)*0.82,TotalPopulation * FractionOfMale*0.82};
    
    /* logMacro(Rprintf("arrayInitialValueForGeneralPopulation_OOO:\n"), 0);
    logDoubleArrayContents(Dimensions_Sex * Dimensions_Age, arrayInitialValueForGeneralPopulation_OOO, 100, 0);*/

  	initializeStateVariablesAtIndexToSubstateSpecificValues(GeneralPopulation_StateIndex, 
                                                            Dimensions_Sex * Dimensions_Age, 
                                                            arrayInitialValueForGeneralPopulation, 
                                                            countStateVariables, 
                                                            stateVariablesForParticleSampled); 
                                                        

    double Initial_RemandPopulationFemale = runif(initial_remandPopulationFemale_min, initial_remandPopulationFemale_max);
    double Initial_SentencedPopulationFemale = runif(initial_SentencedPopulationFemale_min, initial_SentencedPopulationFemale_max);

    double Initial_RemandPopulationMale = runif(initial_remandPopulationMale_min, initial_remandPopulationMale_max);
    double Initial_SentencedPopulationMale = runif(initial_SentencedPopulationMale_min, initial_SentencedPopulationMale_max);

    //double arrayInitialValueForRemandand[Dimensions_Sex] = {Initial_RemandPopulationFemale, Initial_RemandPopulationMale};
    double arrayInitialValueForRemandand[Dimensions_Sex * Dimensions_Age] = {Initial_RemandPopulationFemale*0.12,Initial_RemandPopulationMale*0.12, Initial_RemandPopulationFemale*0.06,Initial_RemandPopulationMale*0.06, Initial_RemandPopulationFemale*0.82, Initial_RemandPopulationMale*0.82};
    //double arrayInitialValueForSentenced[Dimensions_Sex] = {Initial_SentencedPopulationFemale, Initial_SentencedPopulationMale};
    double arrayInitialValueForSentenced[Dimensions_Sex * Dimensions_Age] = {Initial_SentencedPopulationFemale*0.12, Initial_SentencedPopulationMale*0.12, Initial_SentencedPopulationFemale*0.06,Initial_SentencedPopulationMale*0.06,Initial_SentencedPopulationFemale*0.82,Initial_SentencedPopulationMale*0.82};

    initializeStateVariablesAtIndexToSubstateSpecificValues(RemandPopulation_StateIndex, 
                                                            Dimensions_Sex * Dimensions_Age, 
                                                            arrayInitialValueForRemandand, 
                                                            countStateVariables, 
                                                            stateVariablesForParticleSampled); 

    initializeStateVariablesAtIndexToSubstateSpecificValues(SentencedCustody_StateIndex, 
                                                            Dimensions_Sex * Dimensions_Age,
                                                            arrayInitialValueForSentenced, 
                                                            countStateVariables, 
                                                            stateVariablesForParticleSampled);     

    initializeStateVariablesAtIndexToUniformValue(InitialDetentionPopulation_StateIndex, Dimensions_Sex * Dimensions_Age, 10.0, countStateVariables, stateVariablesForParticleSampled);                                                                                                              
    initializeStateVariablesAtIndexToUniformValue(OnBail_StateIndex, Dimensions_Sex * Dimensions_Age, 10.0, countStateVariables, stateVariablesForParticleSampled);
    initializeStateVariablesAtIndexToUniformValue(OnProbation_StateIndex, Dimensions_Sex * Dimensions_Age, 10.0, countStateVariables, stateVariablesForParticleSampled);

    double ArrestHazardRate_PF = runif(arrestHazardRate_PFMin, arrestHazardRate_PFMax);
    stateVariablesForParticleSampled[log_arrestHazard_PF_SEIRStateIndex] = log(ArrestHazardRate_PF);

    double fractionOfInitallyDetainedEnteringRemand_PF = runif(fractionOfInitallyDetainedEnteringRemand_PFMin, fractionOfInitallyDetainedEnteringRemand_PFMax);
    stateVariablesForParticleSampled[logit_fractionOfInitallyDetainedEnteringRemand_PF_StateIndex] = log(fractionOfInitallyDetainedEnteringRemand_PF/(1-fractionOfInitallyDetainedEnteringRemand_PF));;


}

// Integrate the Opioid system from timepoint iTimepoint-1 to iTimepoint 
void integrateODEFromPreviousToCurrentTimepointModelSpecific(int iTimepoint, 
							       int countParameters, 
							       double parameters[countParameters],
							       int countStateVariables, 
							       int countObservationTimepoints,
                                   double arrayTimesForTimepoints[countObservationTimepoints],
                                   // the below is both an INPUT and an OUTPUT variable!
							       double stateVariablesForTimepoint[countObservationTimepoints][countStateVariables])
{

     integrateODEFromPreviousToCurrentTimepointModelSpecificFullyQualified(iTimepoint, 
                                                                           countParameters, 
                                                                           parameters,
                                                                           countStateVariables, 
                                                                           countObservationTimepoints,
                                                                           arrayTimesForTimepoints,
                                                                           integrationTimestep,
                                                                           // the below is both an INPUT and an OUTPUT variable!
                                                                           stateVariablesForTimepoint);
                                                                            
}


void assignStatesToCurrentValuesWithSubscreption(int iStartingIndexToAssign,
                                     int countValuesOfStockToAssign, // total number of stocks with one stock without subscripted
                                     int iTimepointPrevious,
                                     int countStateVariables,    // total stocks in the model
                                     double stateVariablesForTimepoint[iTimepointPrevious][countStateVariables],
                                     double currentStock[countValuesOfStockToAssign])
{
	for (int iValueToAssign = 0; iValueToAssign < countValuesOfStockToAssign; iValueToAssign++)
    	currentStock[iValueToAssign] = stateVariablesForTimepoint[iTimepointPrevious][iStartingIndexToAssign + iValueToAssign];
}

void assignCurrentValuesToStatesWithSubscreption(int iStartingIndexToAssign,
                                     int countValuesOfStockToAssign, // total number of stocks with one stock without subscripted
                                     int iTimepoint,
                                     int countStateVariables,    // total stocks in the model
                                     double stateVariablesForTimepoint[iTimepoint][countStateVariables],
                                     double currentStock[countValuesOfStockToAssign])
{
	for (int iValueToAssign = 0; iValueToAssign < countValuesOfStockToAssign; iValueToAssign++)
    	stateVariablesForTimepoint[iTimepoint][iStartingIndexToAssign + iValueToAssign] = currentStock[iValueToAssign];
}

void calculateFlowWithRate(int countValuesOfStock, // total number of stocks with one stock without subscripted
                   double rateparameter,    // total stocks in the model
                   double currentStockFrom[countValuesOfStock],
                   double currentFlow[countValuesOfStock])
{
	for (int iSubscript = 0; iSubscript < countValuesOfStock; iSubscript++)
    	currentFlow[iSubscript] = currentStockFrom[iSubscript] * rateparameter;
}


void integrateODEFromPreviousToCurrentTimepointModelSpecificFullyQualified(int iTimepoint, 
							       int countParameters, 
							       double parameters[countParameters],
							       int countStateVariables, 
							       int countObservationTimepoints,
                      double arrayTimesForTimepoints[countObservationTimepoints],
                      double integrationTimestepParm,
                      // the below is both an INPUT and an OUTPUT variable!
					  double stateVariablesForTimepoint[countObservationTimepoints][countStateVariables])

{
  double TotalPopulationSize_child = 0;
  double TotalPopulationSize_teen = 0;
  double TotalPopulationSize_old = 0;  
  assert(iTimepoint >= 1); // we need to have a previous timepoint from which to integrate!
  
  int iTimepointPrevious = iTimepoint - 1;
  
  double previousObservationTime = arrayTimesForTimepoints[iTimepointPrevious];             // the starting time of the integration
  double nextObservationTime =  arrayTimesForTimepoints[iTimepoint];                  // the finish time of the integration
  assert(previousObservationTime <= nextObservationTime);                             //ensure that the times associated with the timepoints are ordered
  
  double *stateVariablesForParticleInitial = stateVariablesForTimepoint[iTimepointPrevious];

  double current_GeneralPopulation[Dimensions_Sex * Dimensions_Age];
  double current_InitialDetentionPopulation[Dimensions_Sex * Dimensions_Age];
  double current_RemandPopulation[Dimensions_Sex * Dimensions_Age];
  double current_SentencedCustody[Dimensions_Sex * Dimensions_Age];
  double current_OnBail[Dimensions_Sex * Dimensions_Age];
  double current_OnProbation[Dimensions_Sex * Dimensions_Age];

  //assign the current value of stocks
  assignStatesToCurrentValuesWithSubscreption(GeneralPopulation_StateIndex, Dimensions_Sex * Dimensions_Age, iTimepointPrevious,countStateVariables, stateVariablesForTimepoint, current_GeneralPopulation);
  assignStatesToCurrentValuesWithSubscreption(InitialDetentionPopulation_StateIndex, Dimensions_Sex * Dimensions_Age, iTimepointPrevious,countStateVariables, stateVariablesForTimepoint, current_InitialDetentionPopulation);
  assignStatesToCurrentValuesWithSubscreption(RemandPopulation_StateIndex, Dimensions_Sex * Dimensions_Age, iTimepointPrevious,countStateVariables, stateVariablesForTimepoint, current_RemandPopulation);
  assignStatesToCurrentValuesWithSubscreption(SentencedCustody_StateIndex, Dimensions_Sex * Dimensions_Age, iTimepointPrevious,countStateVariables, stateVariablesForTimepoint, current_SentencedCustody);
  assignStatesToCurrentValuesWithSubscreption(OnBail_StateIndex, Dimensions_Sex * Dimensions_Age, iTimepointPrevious,countStateVariables, stateVariablesForTimepoint, current_OnBail);
  assignStatesToCurrentValuesWithSubscreption(OnProbation_StateIndex, Dimensions_Sex * Dimensions_Age, iTimepointPrevious,countStateVariables, stateVariablesForTimepoint, current_OnProbation);

  double current_log_arrestHazard_PF = stateVariablesForTimepoint[iTimepointPrevious][log_arrestHazard_PF_SEIRStateIndex];
  double current_logit_fractionOfInitallyDetainedEnteringRemand_PF = stateVariablesForTimepoint[iTimepointPrevious][logit_fractionOfInitallyDetainedEnteringRemand_PF_StateIndex];
  


  for (double time = previousObservationTime; time < nextObservationTime; time += integrationTimestep)         // never goes beyond the final time
  {
      //==================================================================================================================================================================================================
      // This is THE KEY PLACE where the state update occurs for integration of the differential equation from the previous observation timepoint to the current observation
      //==================================================================================================================================================================================================

      double dt = min(nextObservationTime - time, integrationTimestep);  // takes into account the amount of time actually remaining; once finished, the next loop to the top will terminate the loop
      
      //doesn't need to add age for this
      double ArrestHazard_PF_male = exp(current_log_arrestHazard_PF);  
      double ArrestHazard_PF_female = ArrestHazard_PF_male * arrestHazardRate_fractionOfFemaleToMale; 
      
      double fractionOfInitallyDetainedEnteringRemand = exp(current_logit_fractionOfInitallyDetainedEnteringRemand_PF)/(1+exp(current_logit_fractionOfInitallyDetainedEnteringRemand_PF));

      // MCMC parameters
      assert(countParameters == 1);
      double logit_fractionOfInitiallyDetainedReturningImmediatelyToCommunity = parameters[logit_fractionOfInitiallyDetainedReturningImmediatelyToCommunity_parameterIndex];
      double fractionOfInitiallyDetainedReturningImmediatelyToCommunity = exp(logit_fractionOfInitiallyDetainedReturningImmediatelyToCommunity)/(1+exp(logit_fractionOfInitiallyDetainedReturningImmediatelyToCommunity));
      /* daily unit
      double meanDurationOfInitialDetention = 0.03 * Month;
      //double arrestHazard = 0.005 / Month; // need double check
      double fractionOfInitallyDetainedEnteringRemand = 0.1;
      //double fractionOfInitiallyDetainedReturningImmediatelyToCommunity = 0.5;
      double fractionSentencedFromRemand = 0.1;
      double courtProcessingPerDay = 40.0 / Month;
      double fractionSentencedFromBail = 0.1;
      double bailToRemandHazard = 0.05 / Month;
      double violateConditionProb = 0.1 ; // in anylogic, it is dynamics variable, but should be parameter.
      double average_days_in_probation = 180.0;
      double average_days_in_remand = 15.0;
      double meanLengthOfSentenceInDays = 36 * Month;
*/   // monthly unit
      double meanDurationOfInitialDetention = 0.03;
      //double arrestHazard = 0.005 / Month; // need double check
      /////double fractionOfInitallyDetainedEnteringRemand = 0.4; // 0.1
      //double fractionOfInitiallyDetainedReturningImmediatelyToCommunity = 0.5;
      double fractionSentencedFromRemand = 0.1;

      double courtProcessingPerDay_male = 250.0; //perMonth 40
      double courtProcessingPerDay_female = 30.0; //perMonth 40

      double fractionSentencedFromBail = 0.1; //0.1
      double bailToRemandHazard = 0.05;
      double violateConditionProb = 0.1 ; // in anylogic, it is dynamics variable, but should be parameter.
      double average_days_in_probation = 180.0/Month; //305
      double average_days_in_remand = 15.0/Month;
      double meanLengthOfSentenceInDays = 36;

      double dischargeWithSupervisionOrRemandFromInitialDetained[Dimensions_Sex * Dimensions_Age];
      double fractionOnRemandVsBail[Dimensions_Sex * Dimensions_Age];
      double fractionAwaitingConclusionWhileInRemand[Dimensions_Sex * Dimensions_Age];
      double totalAwaitingSentencing[Dimensions_Sex * Dimensions_Age];
      double meanTimeUntilSentenced[Dimensions_Sex * Dimensions_Age];
      double courtProcessingPerDayForThoseOnBail[Dimensions_Sex * Dimensions_Age];
      double courtProcessingPerDayForThoseInRemand[Dimensions_Sex * Dimensions_Age];
  


      // calculate dynamic variables
      dischargeWithSupervisionOrRemandFromInitialDetained[Female_child_index] = current_InitialDetentionPopulation[Female_child_index] / meanDurationOfInitialDetention * (1-fractionOfInitiallyDetainedReturningImmediatelyToCommunity);
      dischargeWithSupervisionOrRemandFromInitialDetained[Male_child_index] = (current_InitialDetentionPopulation[Male_child_index] / meanDurationOfInitialDetention)*(1-fractionOfInitiallyDetainedReturningImmediatelyToCommunity);
      dischargeWithSupervisionOrRemandFromInitialDetained[Female_teen_index] = current_InitialDetentionPopulation[Female_teen_index] / meanDurationOfInitialDetention * (1-fractionOfInitiallyDetainedReturningImmediatelyToCommunity);
      dischargeWithSupervisionOrRemandFromInitialDetained[Male_teen_index] = (current_InitialDetentionPopulation[Male_teen_index] / meanDurationOfInitialDetention)*(1-fractionOfInitiallyDetainedReturningImmediatelyToCommunity);
      dischargeWithSupervisionOrRemandFromInitialDetained[Female_old_index] = current_InitialDetentionPopulation[Female_old_index] / meanDurationOfInitialDetention * (1-fractionOfInitiallyDetainedReturningImmediatelyToCommunity);
      dischargeWithSupervisionOrRemandFromInitialDetained[Male_old_index] = (current_InitialDetentionPopulation[Male_old_index] / meanDurationOfInitialDetention)*(1-fractionOfInitiallyDetainedReturningImmediatelyToCommunity);
      
      
      totalAwaitingSentencing[Female_child_index] = current_RemandPopulation[Female_child_index] + current_OnBail[Female_child_index];
      totalAwaitingSentencing[Male_child_index] = current_RemandPopulation[Male_child_index] + current_OnBail[Male_child_index];
      totalAwaitingSentencing[Female_teen_index] = current_RemandPopulation[Female_teen_index] + current_OnBail[Female_teen_index];
      totalAwaitingSentencing[Male_teen_index] = current_RemandPopulation[Male_teen_index] + current_OnBail[Male_teen_index];
      totalAwaitingSentencing[Female_old_index] = current_RemandPopulation[Female_old_index] + current_OnBail[Female_old_index];
      totalAwaitingSentencing[Male_old_index] = current_RemandPopulation[Male_old_index] + current_OnBail[Male_old_index];

      fractionOnRemandVsBail[Female_child_index] = current_RemandPopulation[Female_child_index] / totalAwaitingSentencing[Female_child_index];
      fractionOnRemandVsBail[Male_child_index] = current_RemandPopulation[Male_child_index] / totalAwaitingSentencing[Male_child_index];
      fractionOnRemandVsBail[Female_teen_index] = current_RemandPopulation[Female_teen_index] / totalAwaitingSentencing[Female_teen_index];
      fractionOnRemandVsBail[Male_teen_index] = current_RemandPopulation[Male_teen_index] / totalAwaitingSentencing[Male_teen_index];
      fractionOnRemandVsBail[Female_old_index] = current_RemandPopulation[Female_old_index] / totalAwaitingSentencing[Female_old_index];
      fractionOnRemandVsBail[Male_old_index] = current_RemandPopulation[Male_old_index] / totalAwaitingSentencing[Male_old_index];

      fractionAwaitingConclusionWhileInRemand[Female_child_index] = current_RemandPopulation[Female_child_index] / (current_RemandPopulation[Female_child_index] + current_OnBail[Female_child_index]);
      fractionAwaitingConclusionWhileInRemand[Male_child_index] = current_RemandPopulation[Male_child_index] / (current_RemandPopulation[Male_child_index] + current_OnBail[Male_child_index]);
      fractionAwaitingConclusionWhileInRemand[Female_teen_index] = current_RemandPopulation[Female_teen_index] / (current_RemandPopulation[Female_teen_index] + current_OnBail[Female_teen_index]);
      fractionAwaitingConclusionWhileInRemand[Male_teen_index] = current_RemandPopulation[Male_teen_index] / (current_RemandPopulation[Male_teen_index] + current_OnBail[Male_teen_index]);
      fractionAwaitingConclusionWhileInRemand[Female_old_index] = current_RemandPopulation[Female_old_index] / (current_RemandPopulation[Female_old_index] + current_OnBail[Female_old_index]);
      fractionAwaitingConclusionWhileInRemand[Male_old_index] = current_RemandPopulation[Male_old_index] / (current_RemandPopulation[Male_old_index] + current_OnBail[Male_old_index]);
      
      meanTimeUntilSentenced[Female_child_index] = totalAwaitingSentencing[Female_child_index] / courtProcessingPerDay_female;
      meanTimeUntilSentenced[Male_child_index] = totalAwaitingSentencing[Male_child_index] / courtProcessingPerDay_male;
      meanTimeUntilSentenced[Female_teen_index] = totalAwaitingSentencing[Female_teen_index] / courtProcessingPerDay_female;
      meanTimeUntilSentenced[Male_teen_index] = totalAwaitingSentencing[Male_teen_index] / courtProcessingPerDay_male;
      meanTimeUntilSentenced[Female_old_index] = totalAwaitingSentencing[Female_old_index] / courtProcessingPerDay_female;
      meanTimeUntilSentenced[Male_old_index] = totalAwaitingSentencing[Male_old_index] / courtProcessingPerDay_male;

      courtProcessingPerDayForThoseInRemand[Female_child_index] = fractionOnRemandVsBail[Female_child_index] * courtProcessingPerDay_female;
      courtProcessingPerDayForThoseInRemand[Male_child_index] = fractionOnRemandVsBail[Male_child_index] * courtProcessingPerDay_male;      
      courtProcessingPerDayForThoseInRemand[Female_teen_index] = fractionOnRemandVsBail[Female_teen_index] * courtProcessingPerDay_female;
      courtProcessingPerDayForThoseInRemand[Male_teen_index] = fractionOnRemandVsBail[Male_teen_index] * courtProcessingPerDay_male;     
      courtProcessingPerDayForThoseInRemand[Female_old_index] = fractionOnRemandVsBail[Female_old_index] * courtProcessingPerDay_female;
      courtProcessingPerDayForThoseInRemand[Male_old_index] = fractionOnRemandVsBail[Male_old_index] * courtProcessingPerDay_male;     

      courtProcessingPerDayForThoseOnBail[Female_child_index] = courtProcessingPerDay_female - courtProcessingPerDayForThoseInRemand[Female_child_index];
      courtProcessingPerDayForThoseOnBail[Male_child_index] = courtProcessingPerDay_male - courtProcessingPerDayForThoseInRemand[Male_child_index];
      courtProcessingPerDayForThoseOnBail[Female_teen_index] = courtProcessingPerDay_female - courtProcessingPerDayForThoseInRemand[Female_teen_index];
      courtProcessingPerDayForThoseOnBail[Male_teen_index] = courtProcessingPerDay_male - courtProcessingPerDayForThoseInRemand[Male_teen_index];
      courtProcessingPerDayForThoseOnBail[Female_old_index] = courtProcessingPerDay_female - courtProcessingPerDayForThoseInRemand[Female_old_index];
      courtProcessingPerDayForThoseOnBail[Male_old_index] = courtProcessingPerDay_male - courtProcessingPerDayForThoseInRemand[Male_old_index];      

      #ifdef DEBUG
      logMacro(Rprintf("At end (time %g, time point %d):\n", nextObservationTime, iTimepoint), 0);
      logMacro(Rprintf("the parameter dischargeWithSupervisionOrRemandFromInitialDetained:\n"), 0);
      logDoubleArrayContents(Dimensions_Sex * Dimensions_Age, dischargeWithSupervisionOrRemandFromInitialDetained, 100, 0);
      logMacro(Rprintf("the parameter fractionOnRemandVsBail:\n"), 0);
      logDoubleArrayContents(Dimensions_Sex * Dimensions_Age, fractionOnRemandVsBail, 100, 0);
      logMacro(Rprintf("the parameter fractionAwaitingConclusionWhileInRemand:\n"), 0);
      logDoubleArrayContents(Dimensions_Sex * Dimensions_Age, fractionAwaitingConclusionWhileInRemand, 100, 0);
      logMacro(Rprintf("the parameter totalAwaitingSentencing:\n"), 0);
      logDoubleArrayContents(Dimensions_Sex * Dimensions_Age, totalAwaitingSentencing, 100, 0);
      logMacro(Rprintf("the parameter meanTimeUntilSentenced:\n"), 0);
      logDoubleArrayContents(Dimensions_Sex * Dimensions_Age, meanTimeUntilSentenced, 100, 0);
      logMacro(Rprintf("the parameter courtProcessingPerDayForThoseOnBail:\n"), 0);
      logDoubleArrayContents(Dimensions_Sex * Dimensions_Age, courtProcessingPerDayForThoseOnBail, 100, 0);  
      logMacro(Rprintf("the parameter courtProcessingPerDayForThoseInRemand:\n"), 0);
      logDoubleArrayContents(Dimensions_Sex * Dimensions_Age, courtProcessingPerDayForThoseInRemand, 100, 0);     
      #endif

      // define the flows
      double flow_arrested[Dimensions_Sex * Dimensions_Age];
      double flow_dischargeFromInitialDetensionBackToCommunity[Dimensions_Sex * Dimensions_Age];
      double flow_enteringRemand[Dimensions_Sex * Dimensions_Age];
      double flow_becomingSentencedOffender[Dimensions_Sex * Dimensions_Age];
      double flow_returningToCommunity[Dimensions_Sex * Dimensions_Age];
      double flow_sentencedOffendersFinishingSentence[Dimensions_Sex * Dimensions_Age];
      double flow_bailToCommunity[Dimensions_Sex * Dimensions_Age];
      double flow_enteringBail[Dimensions_Sex * Dimensions_Age];
      double flow_violatingConditions[Dimensions_Sex * Dimensions_Age];
      double flow_sentencedAfterBail[Dimensions_Sex * Dimensions_Age];
      double flow_releaseFromProbation[Dimensions_Sex * Dimensions_Age];
      double flow_releaseToProbabtion[Dimensions_Sex * Dimensions_Age];
      double flow_violatedProbationConditions[Dimensions_Sex * Dimensions_Age];
      
      double flow_ageout_child_GeneralPoluation[Dimensions_Sex];
      double flow_ageout_child_InitialDetentionPopulation[Dimensions_Sex];
      double flow_ageout_child_RemandPopulation[Dimensions_Sex];
      double flow_ageout_child_SentencedCustody[Dimensions_Sex];
      double flow_ageout_child_OnBail[Dimensions_Sex];
      double flow_ageout_child_OnProbation[Dimensions_Sex];

      double flow_ageout_teen_GeneralPoluation[Dimensions_Sex];
      double flow_ageout_teen_InitialDetentionPopulation[Dimensions_Sex];
      double flow_ageout_teen_RemandPopulation[Dimensions_Sex];
      double flow_ageout_teen_SentencedCustody[Dimensions_Sex];
      double flow_ageout_teen_OnBail[Dimensions_Sex];
      double flow_ageout_teen_OnProbation[Dimensions_Sex];   

      // calculate the flows
      flow_arrested[Female_child_index] = current_GeneralPopulation[Female_child_index] * ArrestHazard_PF_female;
      flow_arrested[Male_child_index] = current_GeneralPopulation[Male_child_index] * ArrestHazard_PF_male;
      calculateFlowWithRate(Dimensions_Sex * Dimensions_Age, (fractionOfInitiallyDetainedReturningImmediatelyToCommunity / meanDurationOfInitialDetention), current_InitialDetentionPopulation, flow_dischargeFromInitialDetensionBackToCommunity);
      calculateFlowWithRate(Dimensions_Sex * Dimensions_Age, fractionOfInitallyDetainedEnteringRemand, dischargeWithSupervisionOrRemandFromInitialDetained, flow_enteringRemand);
      
      calculateFlowWithRate(Dimensions_Sex * Dimensions_Age, fractionSentencedFromRemand, courtProcessingPerDayForThoseInRemand, flow_becomingSentencedOffender);
      calculateFlowWithRate(Dimensions_Sex * Dimensions_Age, (1.0 - fractionSentencedFromRemand), courtProcessingPerDayForThoseInRemand, flow_returningToCommunity);
      
      for (int iSubscript = 0; iSubscript < Dimensions_Sex; iSubscript ++)
      {
          flow_ageout_child_GeneralPoluation[iSubscript] = current_GeneralPopulation[iSubscript] * aging_out_child;
          flow_ageout_child_InitialDetentionPopulation[iSubscript] = current_InitialDetentionPopulation[iSubscript] * aging_out_child;
          flow_ageout_child_RemandPopulation[iSubscript] = current_RemandPopulation[iSubscript] * aging_out_child;
          flow_ageout_child_SentencedCustody[iSubscript] = current_SentencedCustody[iSubscript] * aging_out_child;
          flow_ageout_child_OnBail[iSubscript] = current_OnBail[iSubscript] * aging_out_child;
          flow_ageout_child_OnProbation[iSubscript] = current_OnProbation[iSubscript] * aging_out_child;

          flow_ageout_teen_GeneralPoluation[iSubscript] = current_GeneralPopulation[iSubscript] * aging_out_teen;
          flow_ageout_teen_InitialDetentionPopulation[iSubscript] = current_InitialDetentionPopulation[iSubscript] * aging_out_teen;
          flow_ageout_teen_RemandPopulation[iSubscript] = current_RemandPopulation[iSubscript] * aging_out_teen;
          flow_ageout_teen_SentencedCustody[iSubscript] = current_SentencedCustody[iSubscript] * aging_out_teen;
          flow_ageout_teen_OnBail[iSubscript] = current_OnBail[iSubscript] * aging_out_teen;
          flow_ageout_teen_OnProbation[iSubscript] = current_OnProbation[iSubscript] * aging_out_teen;
      }

      
      #ifdef  DEBUG_D
      logMacro(Rprintf("ArrestHazard_PF_female=%g\n", ArrestHazard_PF_female), 0);
      logMacro(Rprintf("ArrestHazard_PF_male=%g\n", ArrestHazard_PF_male), 0);
      logMacro(Rprintf("fractionOfInitallyDetainedEnteringRemand=%g\n", fractionOfInitallyDetainedEnteringRemand), 0);

      #endif  


      /*
      flow_becomingSentencedOffender[Female_child_index] = courtProcessingPerDayForThoseInRemand[Female_child_index] * fractionSentencedFromRemand_female;
      flow_returningToCommunity[Female_child_index] = courtProcessingPerDayForThoseInRemand[Female_child_index] * (1.0 - fractionSentencedFromRemand_female);
      flow_becomingSentencedOffender[Male_child_index] = courtProcessingPerDayForThoseInRemand[Male_child_index] * fractionSentencedFromRemand_male;
      flow_returningToCommunity[Male_child_index] = courtProcessingPerDayForThoseInRemand[Male_child_index] * (1.0 - fractionSentencedFromRemand_male);
      */
      calculateFlowWithRate(Dimensions_Sex * Dimensions_Age, (1.0 / meanLengthOfSentenceInDays),  current_SentencedCustody, flow_sentencedOffendersFinishingSentence);
      calculateFlowWithRate(Dimensions_Sex * Dimensions_Age, (1.0 - fractionSentencedFromBail), courtProcessingPerDayForThoseOnBail, flow_bailToCommunity);
      calculateFlowWithRate(Dimensions_Sex * Dimensions_Age, (1.0 - fractionOfInitallyDetainedEnteringRemand), dischargeWithSupervisionOrRemandFromInitialDetained, flow_enteringBail);
      calculateFlowWithRate(Dimensions_Sex * Dimensions_Age, bailToRemandHazard, current_OnBail, flow_violatingConditions);
      calculateFlowWithRate(Dimensions_Sex * Dimensions_Age, fractionSentencedFromBail, courtProcessingPerDayForThoseOnBail, flow_sentencedAfterBail);
      calculateFlowWithRate(Dimensions_Sex * Dimensions_Age, (1.0 / average_days_in_probation), current_OnProbation, flow_releaseFromProbation);
      calculateFlowWithRate(Dimensions_Sex * Dimensions_Age, (1.0 / average_days_in_remand), current_RemandPopulation, flow_releaseToProbabtion);
      calculateFlowWithRate(Dimensions_Sex * Dimensions_Age, (violateConditionProb / Year), current_OnProbation, flow_violatedProbationConditions);
      #ifdef DEBUG_flow
      logMacro(Rprintf("At end (time %g, time point %d):\n", nextObservationTime, iTimepoint), 0);
      logMacro(Rprintf("the flow flow_arrested:\n"), 0);
      logDoubleArrayContents(Dimensions_Sex * Dimensions_Age, flow_arrested, 100, 0);
      logMacro(Rprintf("the flow flow_dischargeFromInitialDetensionBackToCommunity:\n"), 0);
      logDoubleArrayContents(Dimensions_Sex * Dimensions_Age, flow_dischargeFromInitialDetensionBackToCommunity, 100, 0);
      logMacro(Rprintf("the flow flow_enteringRemand:\n"), 0);
      logDoubleArrayContents(Dimensions_Sex * Dimensions_Age, flow_enteringRemand, 100, 0);
      logMacro(Rprintf("the flow flow_becomingSentencedOffender:\n"), 0);
      logDoubleArrayContents(Dimensions_Sex * Dimensions_Age, flow_becomingSentencedOffender, 100, 0);
      logMacro(Rprintf("the flow flow_returningToCommunity:\n"), 0);
      logDoubleArrayContents(Dimensions_Sex * Dimensions_Age, flow_returningToCommunity, 100, 0);
      logMacro(Rprintf("the flow flow_sentencedOffendersFinishingSentence:\n"), 0);
      logDoubleArrayContents(Dimensions_Sex * Dimensions_Age, flow_sentencedOffendersFinishingSentence, 100, 0);  
      logMacro(Rprintf("the flow flow_bailToCommunity:\n"), 0);
      logDoubleArrayContents(Dimensions_Sex * Dimensions_Age, flow_bailToCommunity, 100, 0); 

      logMacro(Rprintf("the flow flow_enteringBail:\n"), 0);
      logDoubleArrayContents(Dimensions_Sex * Dimensions_Age, flow_enteringBail, 100, 0);
      logMacro(Rprintf("the flow flow_violatingConditions:\n"), 0);
      logDoubleArrayContents(Dimensions_Sex * Dimensions_Age, flow_violatingConditions, 100, 0);
      logMacro(Rprintf("the flow flow_sentencedAfterBail:\n"), 0);
      logDoubleArrayContents(Dimensions_Sex * Dimensions_Age, flow_sentencedAfterBail, 100, 0);
      logMacro(Rprintf("the flow flow_releaseFromProbation:\n"), 0);
      logDoubleArrayContents(Dimensions_Sex * Dimensions_Age, flow_releaseFromProbation, 100, 0);
      logMacro(Rprintf("the flow flow_releaseToProbabtion:\n"), 0);
      logDoubleArrayContents(Dimensions_Sex * Dimensions_Age, flow_releaseToProbabtion, 100, 0);
      logMacro(Rprintf("the flow flow_violatedProbationConditions:\n"), 0);
      logDoubleArrayContents(Dimensions_Sex * Dimensions_Age, flow_violatedProbationConditions, 100, 0);  

      #endif

      // update the stocks: child = -aging_out + in - out (sentenced)
      //double TotalPopulationSize_child = 0;
      for (int iSubscript = 0; iSubscript < Dimensions_Sex; iSubscript ++)
      {
      	  current_GeneralPopulation[iSubscript] += (+ flow_returningToCommunity[iSubscript]
      	                                           + flow_releaseFromProbation[iSubscript]
      	                                           + flow_dischargeFromInitialDetensionBackToCommunity[iSubscript]
      	                                           + flow_sentencedOffendersFinishingSentence[iSubscript]
      	                                           + flow_bailToCommunity[iSubscript]
      	                                           - flow_arrested[iSubscript]
                                                   - flow_ageout_child_GeneralPoluation[iSubscript]) * dt;

      	  current_InitialDetentionPopulation[iSubscript] += (+ flow_arrested[iSubscript]
      	                                                    - flow_dischargeFromInitialDetensionBackToCommunity[iSubscript]
      	                                                    - flow_enteringRemand[iSubscript]
      	                                                    - flow_enteringBail[iSubscript]
                                                            - flow_ageout_child_InitialDetentionPopulation[iSubscript]) * dt;

      	  current_RemandPopulation[iSubscript] += (+ flow_enteringRemand[iSubscript]  //here!!! = -aging_out + in - out (sentenced)
      	                                          + flow_violatingConditions[iSubscript]
      	                                          - flow_releaseToProbabtion[iSubscript]
      	                                          - flow_returningToCommunity[iSubscript]
      	                                          - flow_becomingSentencedOffender[iSubscript]
                                                  - flow_ageout_child_RemandPopulation[iSubscript]) * dt;

      	  current_SentencedCustody[iSubscript] += (+ flow_violatedProbationConditions[iSubscript]
      	                                          + flow_becomingSentencedOffender[iSubscript]
      	                                          + flow_sentencedAfterBail[iSubscript]
      	                                          - flow_sentencedOffendersFinishingSentence[iSubscript]
                                                  - flow_ageout_child_SentencedCustody[iSubscript]) * dt;

      	  current_OnBail[iSubscript] += (+ flow_enteringBail[iSubscript]
      	                                - flow_bailToCommunity[iSubscript]
      	                                - flow_violatingConditions[iSubscript]
      	                                - flow_sentencedAfterBail[iSubscript]
                                        - flow_ageout_child_OnBail[iSubscript]) * dt;

      	  current_OnProbation[iSubscript] += (+ flow_releaseToProbabtion[iSubscript]
      	                                     - flow_violatedProbationConditions[iSubscript]
      	                                     - flow_releaseFromProbation[iSubscript]
                                             - flow_ageout_child_OnProbation[iSubscript]) * dt;
          TotalPopulationSize_child += current_GeneralPopulation[iSubscript] 
                                 + current_InitialDetentionPopulation[iSubscript] 
                                 + current_RemandPopulation[iSubscript] 
                                 + current_OnProbation[iSubscript] 
                                 + current_OnBail[iSubscript]
                                 + current_SentencedCustody[iSubscript];

      }




      // update the stocks: teen = + aging_in - aging_out + in - out (sentenced)
      for (int iSubscript = Dimensions_Sex; iSubscript < 2*Dimensions_Sex; iSubscript ++)
      {
      	  current_GeneralPopulation[iSubscript] += (+ flow_returningToCommunity[iSubscript]
      	                                           + flow_releaseFromProbation[iSubscript]
      	                                           + flow_dischargeFromInitialDetensionBackToCommunity[iSubscript]
      	                                           + flow_sentencedOffendersFinishingSentence[iSubscript]
      	                                           + flow_bailToCommunity[iSubscript]
      	                                           - flow_arrested[iSubscript]
                                                   + flow_ageout_child_GeneralPoluation[iSubscript]
                                                   - flow_ageout_teen_GeneralPoluation[iSubscript]) * dt;

      	  current_InitialDetentionPopulation[iSubscript] += (+ flow_arrested[iSubscript]
      	                                                    - flow_dischargeFromInitialDetensionBackToCommunity[iSubscript]
      	                                                    - flow_enteringRemand[iSubscript]
      	                                                    - flow_enteringBail[iSubscript]
                                                            + flow_ageout_child_InitialDetentionPopulation[iSubscript]
                                                            - flow_ageout_teen_InitialDetentionPopulation[iSubscript]) * dt;


      	  current_RemandPopulation[iSubscript] += (+ flow_enteringRemand[iSubscript]
      	                                          + flow_violatingConditions[iSubscript]
      	                                          - flow_releaseToProbabtion[iSubscript]
      	                                          - flow_returningToCommunity[iSubscript]
      	                                          - flow_becomingSentencedOffender[iSubscript]
                                                  + flow_ageout_child_RemandPopulation[iSubscript]
                                                  - flow_ageout_teen_RemandPopulation[iSubscript]) * dt;

      	  current_SentencedCustody[iSubscript] += (+ flow_violatedProbationConditions[iSubscript]
      	                                          + flow_becomingSentencedOffender[iSubscript]
      	                                          + flow_sentencedAfterBail[iSubscript]
      	                                          - flow_sentencedOffendersFinishingSentence[iSubscript]
                                                  + flow_ageout_child_SentencedCustody[iSubscript]
                                                  - flow_ageout_teen_SentencedCustody[iSubscript]) * dt;
      	  current_OnBail[iSubscript] += (+ flow_enteringBail[iSubscript]
      	                                - flow_bailToCommunity[iSubscript]
      	                                - flow_violatingConditions[iSubscript]
      	                                - flow_sentencedAfterBail[iSubscript]
                                        + flow_ageout_child_OnBail[iSubscript]
                                        - flow_ageout_teen_OnBail[iSubscript]) * dt;

      	  current_OnProbation[iSubscript] += (+ flow_releaseToProbabtion[iSubscript]
      	                                     - flow_violatedProbationConditions[iSubscript]
      	                                     - flow_releaseFromProbation[iSubscript]
                                             + flow_ageout_child_OnProbation[iSubscript]
                                             - flow_ageout_teen_OnProbation[iSubscript]) * dt;


      }

      // update the stocks: old = + aging_in + in - out (sentenced)
      for (int iSubscript = 2*Dimensions_Sex; iSubscript < 3*Dimensions_Sex; iSubscript ++)
      {
      	  current_GeneralPopulation[iSubscript] += (+ flow_returningToCommunity[iSubscript]
      	                                           + flow_releaseFromProbation[iSubscript]
      	                                           + flow_dischargeFromInitialDetensionBackToCommunity[iSubscript]
      	                                           + flow_sentencedOffendersFinishingSentence[iSubscript]
      	                                           + flow_bailToCommunity[iSubscript]
      	                                           - flow_arrested[iSubscript]
                                                   + flow_ageout_teen_GeneralPoluation[iSubscript]) * dt;

      	  current_InitialDetentionPopulation[iSubscript] += (+ flow_arrested[iSubscript]
      	                                                    - flow_dischargeFromInitialDetensionBackToCommunity[iSubscript]
      	                                                    - flow_enteringRemand[iSubscript]
      	                                                    - flow_enteringBail[iSubscript]
                                                            + flow_ageout_teen_InitialDetentionPopulation[iSubscript]) * dt;

      	  current_RemandPopulation[iSubscript] += (+ flow_enteringRemand[iSubscript]
      	                                          + flow_violatingConditions[iSubscript]
      	                                          - flow_releaseToProbabtion[iSubscript]
      	                                          - flow_returningToCommunity[iSubscript]
      	                                          - flow_becomingSentencedOffender[iSubscript]
                                                  + flow_ageout_teen_RemandPopulation[iSubscript]) * dt;

      	  current_SentencedCustody[iSubscript] += (+ flow_violatedProbationConditions[iSubscript]
      	                                          + flow_becomingSentencedOffender[iSubscript]
      	                                          + flow_sentencedAfterBail[iSubscript]
      	                                          - flow_sentencedOffendersFinishingSentence[iSubscript]
                                                  + flow_ageout_teen_SentencedCustody[iSubscript]) * dt;

      	  current_OnBail[iSubscript] += (+ flow_enteringBail[iSubscript]
      	                                - flow_bailToCommunity[iSubscript]
      	                                - flow_violatingConditions[iSubscript]
      	                                - flow_sentencedAfterBail[iSubscript]
                                        + flow_ageout_teen_OnBail[iSubscript]) * dt;

      	  current_OnProbation[iSubscript] += (+ flow_releaseToProbabtion[iSubscript]
      	                                     - flow_violatedProbationConditions[iSubscript]
      	                                     - flow_releaseFromProbation[iSubscript]
                                             + flow_ageout_teen_OnProbation[iSubscript]) * dt;


      }
    



      double perturb_arrestHazard_PF = rnorm(0, arrestHazard_RandomWalkStdDev);
      current_log_arrestHazard_PF += (perturb_arrestHazard_PF)*dt; 

      double perturb_fractionOfInitallyDetainedEnteringRemand_PF = rnorm(0, fractionOfInitallyDetainedEnteringRemand_RandomWalkStdDev);
      current_logit_fractionOfInitallyDetainedEnteringRemand_PF += (perturb_fractionOfInitallyDetainedEnteringRemand_PF)*dt; 


    } 

      assignCurrentValuesToStatesWithSubscreption(GeneralPopulation_StateIndex, Dimensions_Sex * Dimensions_Age, iTimepoint,countStateVariables, stateVariablesForTimepoint, current_GeneralPopulation);
      assignCurrentValuesToStatesWithSubscreption(InitialDetentionPopulation_StateIndex, Dimensions_Sex * Dimensions_Age, iTimepoint,countStateVariables, stateVariablesForTimepoint, current_InitialDetentionPopulation);
      assignCurrentValuesToStatesWithSubscreption(RemandPopulation_StateIndex, Dimensions_Sex * Dimensions_Age, iTimepoint,countStateVariables, stateVariablesForTimepoint, current_RemandPopulation);
      assignCurrentValuesToStatesWithSubscreption(SentencedCustody_StateIndex, Dimensions_Sex * Dimensions_Age, iTimepoint,countStateVariables, stateVariablesForTimepoint, current_SentencedCustody);
      assignCurrentValuesToStatesWithSubscreption(OnBail_StateIndex, Dimensions_Sex * Dimensions_Age, iTimepoint,countStateVariables, stateVariablesForTimepoint, current_OnBail);
      assignCurrentValuesToStatesWithSubscreption(OnProbation_StateIndex, Dimensions_Sex * Dimensions_Age, iTimepoint,countStateVariables, stateVariablesForTimepoint, current_OnProbation);
      //assignCurrentValuesToStatesWithSubscreption(OnProbation_StateIndex, Dimensions_Sex, iTimepoint,countStateVariables, stateVariablesForTimepoint, current_OnProbation);

      stateVariablesForTimepoint[iTimepoint][log_arrestHazard_PF_SEIRStateIndex] = current_log_arrestHazard_PF;
      stateVariablesForTimepoint[iTimepoint][logit_fractionOfInitallyDetainedEnteringRemand_PF_StateIndex] = current_logit_fractionOfInitallyDetainedEnteringRemand_PF;

      double TotalPopulationSize = 0;
      for(int iSubscript = 0; iSubscript < Dimensions_Sex * Dimensions_Age; iSubscript++)
      {
          TotalPopulationSize += current_GeneralPopulation[iSubscript] 
                                 + current_InitialDetentionPopulation[iSubscript] 
                                 + current_RemandPopulation[iSubscript] 
                                 + current_OnProbation[iSubscript] 
                                 + current_OnBail[iSubscript]
                                 + current_SentencedCustody[iSubscript];
      }

      #ifdef DEBUG_stock
      logMacro(Rprintf("At end (time %g, time point %d):\n", nextObservationTime, iTimepoint), 0);
      logMacro(Rprintf("the Stock current_GeneralPopulation:\n"), 0);
      logDoubleArrayContents(Dimensions_Sex, current_GeneralPopulation, 100, 0);
      logMacro(Rprintf("the Stock current_InitialDetentionPopulation:\n"), 0);
      logDoubleArrayContents(Dimensions_Sex, current_InitialDetentionPopulation, 100, 0);
      logMacro(Rprintf("the Stock current_RemandPopulation:\n"), 0);
      logDoubleArrayContents(Dimensions_Sex, current_RemandPopulation, 100, 0);
      logMacro(Rprintf("the Stock current_SentencedCustody:\n"), 0);
      logDoubleArrayContents(Dimensions_Sex, current_SentencedCustody, 100, 0);
      logMacro(Rprintf("the Stock current_OnBail:\n"), 0);
      logDoubleArrayContents(Dimensions_Sex, current_OnBail, 100, 0);
      logMacro(Rprintf("the Stock current_OnProbation:\n"), 0);
      logDoubleArrayContents(Dimensions_Sex, current_OnProbation, 100, 0);
      logMacro(Rprintf("The total Population is %g:\n", TotalPopulationSize), 0);
      logMacro(Rprintf("The total Population_child is %g:\n", TotalPopulationSize_child), 0);      
      #endif
  }


#define rCountOfNonCasesReportedDispersionParameter_female (2) //5
#define rCountOfNonCasesReportedDispersionParameter_male (2) //10

double computeLikelihoodOfObservationGivenModelStateModelSpecific( int iTimePoint,
                  int countObservables, 
								  double observationsForTimepoint[countObservables],
								  int countParameters,
								  double parameters[countParameters],
								  int countStateVariables, 
								  double stateVariablesForParticle[countStateVariables])
{

  assert(countObservables == 4);     // for this particular model, there is one observable
  assert(countStateVariables == 14); 

  double countRemandFromModel_male = stateVariablesForParticle[RemandPopulation_StateIndex + Male_child_index];
  double countRemandFromModel_female = stateVariablesForParticle[RemandPopulation_StateIndex + Female_child_index];
  double countStencedCustodyFromModel_male = stateVariablesForParticle[SentencedCustody_StateIndex + Male_child_index];
  double countStencedCustodyFromModel_female = stateVariablesForParticle[SentencedCustody_StateIndex + Female_child_index];

  double countRemandEmpirical_male = observationsForTimepoint[Empirical_Remand_Male_Index];
  double countRemandEmpirical_female = observationsForTimepoint[Empirical_Remand_Female_Index];
  double countSentencedCustodyEmpirical_male = observationsForTimepoint[Empirical_SentencedCustody_Male_Index];
  double countSentencedCustodyEmpirical_female = observationsForTimepoint[Empirical_SentencedCustody_Female_Index];

  double probFailure_Remand_male = rCountOfNonCasesReportedDispersionParameter_male / (countRemandFromModel_male + rCountOfNonCasesReportedDispersionParameter_male);
  double probFailure_Remand_female = rCountOfNonCasesReportedDispersionParameter_female / (countRemandFromModel_female + rCountOfNonCasesReportedDispersionParameter_female);
  double probFailure_SentencedCustody_male = rCountOfNonCasesReportedDispersionParameter_male / (countStencedCustodyFromModel_male + rCountOfNonCasesReportedDispersionParameter_male);
  double probFailure_SentencedCustody_female = rCountOfNonCasesReportedDispersionParameter_female / (countStencedCustodyFromModel_female + rCountOfNonCasesReportedDispersionParameter_female);

  int isLog = 0;
  
  #ifdef  DEBUG_LIKE
  //****TODO AFTER DEUBGGING:  change logging level below back to 5
  logAtLoggingRankInternal("#10\n", 0);
  logMacro(Rprintf("countRemandEmpirical_male=%g;countRemandFromModel_male=%g; probFailure_Remand_male=%g\n", countRemandEmpirical_male, countRemandFromModel_male, probFailure_Remand_male), 0);
  logMacro(Rprintf("countRemandEmpirical_female=%g;countRemandFromModel_female=%g; probFailure_Remand_female=%g\n", countRemandEmpirical_female, countRemandFromModel_female, probFailure_Remand_female), 0);
  logMacro(Rprintf("countSentencedCustodyEmpirical_male=%g;countStencedCustodyFromModel_male=%g; probFailure_SentencedCustody_male=%g\n", countSentencedCustodyEmpirical_male, countStencedCustodyFromModel_male, probFailure_SentencedCustody_male), 0);
  logMacro(Rprintf("countSentencedCustodyEmpirical_female=%g;countStencedCustodyFromModel_female=%g; probFailure_SentencedCustody_female=%g\n", countSentencedCustodyEmpirical_female, countStencedCustodyFromModel_female, probFailure_SentencedCustody_female), 0);

  #endif  

  double likelihood_remand_male = dnbinom(countRemandEmpirical_male, rCountOfNonCasesReportedDispersionParameter_male, probFailure_Remand_male, isLog); 
  double likelihood_remand_female = dnbinom(countRemandEmpirical_female, rCountOfNonCasesReportedDispersionParameter_female, probFailure_Remand_female, isLog); 
  double likelihood_sentencedCustody_male = dnbinom(countSentencedCustodyEmpirical_male, rCountOfNonCasesReportedDispersionParameter_male, probFailure_SentencedCustody_male, isLog); 
  double likelihood_sentencedCustody_female = dnbinom(countSentencedCustodyEmpirical_female, rCountOfNonCasesReportedDispersionParameter_female, probFailure_SentencedCustody_female, isLog); 
  logMacro(Rprintf("likelihood_remand_male %d,",likelihood_remand_male), 0);//48206688
  //logMacro(Rprintf("rCountOfNonCasesReportedDispersionParameter_male %d,",rCountOfNonCasesReportedDispersionParameter_male), 0);//2

//  logMacro(Rprintf("probFailure_Remand_male %d,",probFailure_Remand_male), 0);  //770029408

  double likelihood = likelihood_remand_female * likelihood_remand_male * likelihood_sentencedCustody_female * likelihood_sentencedCustody_male  * pow(10, 10) ;

  #ifdef  DEBUG_LIKELIHOOD
  logMacro(Rprintf("likelihood_remand_male=%g\n", likelihood_remand_male), 0);
  logMacro(Rprintf("likelihood_remand_female=%g\n", likelihood_remand_female), 0);
  logMacro(Rprintf("likelihood_sentencedCustody_female=%g\n", likelihood_sentencedCustody_female), 0);
  logMacro(Rprintf("likelihood_sentencedCustody_male=%g\n", likelihood_sentencedCustody_male), 0);
  logMacro(Rprintf("likelihood=%g\n", likelihood), 0666);
  #endif  

  return(likelihood);

}

double computePriorModelSpecific(int countParameters, 
                                 double parameters[countParameters])
{
    return(1.0L);    // using uninformative prior 
}



