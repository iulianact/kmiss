Simulation code for manuscript:
Causal inference in survival analysis under deterministic missingness of confounders in register data

Note: correspondence between code and manuscript names and abbreviations:
Scenarios:
	Scenario 1 in the code corresponds to Scenario I in the manuscript.
	Scenario 3 in the code corresponds to Scenario II in the manuscript.
Imputation models:
	mice3 in the code corresponds to mice4 in the manuscript.
	mice4 in the code corresponds to mice3 in the manuscript.

Scripts description:
- functions.R: define auxiliary functions
- set_params.R: set general simulation parameters
- sc1.R, sc3.R: set scenario-specific simulation parameters 
- simulation_gen.R: generate full data and imputed data
- simulation_fbc_std.R: full, ignore and compl analysis with standardisation
- simulation_fbc_PS.R: full, ignore and compl analysis with propensity score weighting
- simulation_mice_std.R: imputation with standardisation
- simulation_mice_PS.R: imputation with propensity score weighting
- simulation_mice3_other.R: imputation with other methods.