from run_model import run_model
from joblib import Parallel, delayed

# Define the inputs for the function
inputs = [7]
outcome_type = "births" 
cats = ['total']
dists = ['NB'] # Poisson or NB
missing_flags = [True]
# disp_params = [1e-4, 1e-3]
disp_params = [1e-4]
## placebo_times = ["2020-05-01"]
placebo_times = ["2016-03-01", "2016-05-01", "2016-07-01", "2016-09-01", "2016-11-01",
                 "2017-03-01", "2017-05-01", "2017-07-01", "2017-09-01", "2017-11-01",
                 "2018-03-01", "2018-05-01", "2018-07-01", "2018-09-01", "2018-11-01",
                 "2019-03-01", "2019-05-01", "2019-07-01", "2019-09-01", "2019-11-01",
                 "2020-03-01"]
placebo_states = [None]
sample_disp = False

args = [(dist, cat, rank, m, disp, p, tm) for dist in dists for rank in inputs for cat in cats 
        for m in missing_flags for disp in disp_params for p in placebo_states 
        for tm in placebo_times]
# Run the function in parallel
results = Parallel(n_jobs=100)(delayed(run_model)(dist=i[0], outcome_type=outcome_type, cat_name=i[1], 
                                                  rank=i[2], missingness=i[3], 
                                                  disp_param=i[4], sample_disp=sample_disp, placebo_state=i[5], 
                                                  placebo_time = i[6], results_file_suffix="placebo_" + i[6], 
                                                  num_chains=4, num_samples=2500, 
                                                  num_warmup=1000, thinning=10) for i in args)

for date in placebo_times:
    print(date)
    run_model(dist = "NB", outcome_type="births", cat_name="total", rank=7, 
            missingness=True, disp_param=1e-4, sample_disp=False, 
            placebo_state=None, placebo_time=date, model_treated=True, 
            results_file_suffix="tmp", num_chains=1, num_samples=100, num_warmup=100, thinning=1)