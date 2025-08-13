from run_model import run_model
from clean_monthly_birth_data import subgroup_definitions
from joblib import Parallel, delayed

# Define the inputs for the function
inputs = [6, 7, 8, 9, 10, 11, 12]
cats = list(subgroup_definitions[outcome_type].keys())
outcome_type = "births" 
dists = ['NB'] # Poisson or NB
missing_flags = [True]
# disp_params = [1e-4, 1e-3]
disp_params = [1e-4]
## placebo_times = ["2020-05-01"]
placebo_times = [None]
placebo_states = [None]
sample_disp = False

args = [(dist, cat, rank, m, disp, p, tm) for dist in dists for rank in inputs for cat in cats 
        for m in missing_flags for disp in disp_params for p in placebo_states 
        for tm in placebo_times]
# Run the function in parallel
results = Parallel(n_jobs=100)(delayed(run_model)(dist=i[0], outcome_type=outcome_type, cat_name=i[1], 
                                                  rank=i[2], missingness=i[3], 
                                                  disp_param=i[4], sample_disp=sample_disp, placebo_state=i[5], 
                                                  placebo_time = i[6], results_file_suffix="sensitivity", 
                                                  num_chains=4, num_samples=2500, 
                                                  dobbs_donor_sensitivity=True,
                                                  num_warmup=1000, thinning=10) for i in args)

