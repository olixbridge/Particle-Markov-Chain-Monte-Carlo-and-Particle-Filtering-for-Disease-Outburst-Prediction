| What Olivia did |   |
| --- | --- |
| **AGEGROUP** Added AgeGroup structure to the original model     |
-
  - child(0-12), teen(12-18), old(18+)
  - From 14 stocks to 38 stocks
  - From 2 dimension (gender) to 6 dimension (gender\*age)
  - Modified remandModel.c, remandModel.h files
  - Added &quot;ageout\_child&quot; and &quot;ageout\_teen&quot; flow and modified the integration from 1 section to 3 section (3 forloops).
 |
| Updated likelihood function to log likelihood function |
-
  - Solved the &quot;newWeight = nan&quot; problem, it was because one of the likelihood = nan.
 |
| Updated graphs in R |
-
  - Changed the graphs from female, male to female\_child, female\_teen, female\_old, male\_child, male\_teen, male\_old
  - Modified RInterfaceRemand file
  - Note: Indices are off in the R file: the index in R should be index in remandModel.h+1
 |
|   |   |
| **VARY** |
-
  - arrest rate varies between a minimum (0.015) and a maximum value (0.03)
 |
|   |   |



| **Future work** |   |
| --- | --- |
| Add new empirical data | Column1:year (for double check)Column2:inPolice teen (detention)Column3:inPolice old (detention) Remand and sentenced custody for child = 0 |
| Double check &quot;countTimePoints&quot;\&lt;- 133 | I changed countTimePoints from 299 to 133 because: stopifnot(length(arrayTimesForTimepoints) \* countObservables == length(arrayObservedValuesYByTimepointAndObservable)) |
| Adjust initial value for stocks | Initial detention value |
| How to run PMCMC  | **Run PMCMC:** In terminal: Set the path to be in the foldercd + path docker-compose up -d docker exec -it rstudio bash cd code/ ulimit -s unlimited Rscript RInterfaceRemand.R |
| How to fix PMCMC | ls -l docker ps docker container rm -f rstudio docker exec -it rstudio bash  **If modified the C file:** make cleanmake  |





| **PF** |   |
| --- | --- |
| IZ\_V15  |
-
  - Has measles, chicken pox and influenza with respective empirical values.
  - Immunity is fixed
  - Added Calibration (Optimization Experiment)
  - I ignored the histograms and plots so that calibration can run faster
 |
|   |   |
|   |   |
| Rubella V1 |
-
  - Has measles, chicken pox and Rubella with respective empirical values.
  - Has calibration and histograms
 |
| Rubella V2 |
-
  - Has calibration and but no histograms
 |
