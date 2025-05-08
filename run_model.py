# import os
# os.chdir("birthrate_mtgp")
from jax import numpy as jnp
import numpy as np
import numpyro.distributions as dist
import jax.numpy as jnp
import numpyro
from numpyro.handlers import scope

from models.panel_nmf_model import model
from models.utils import missingness_adjustment
from numpyro_to_draws_df_csv import dict_to_tidybayes

import pandas as pd

## Some defaults for testing
dist = "NB"
outcome_type = "births"
cat_name = "total"
rank = 5
sample_disp = False
missingness=True
disp_param = 1e-4
model_treated = True
placebo_time = "2019-03-01"
num_chains = 1
num_samples = 100
num_warmup=100
start_time = '2016-01-01'
end_time = '2024-01-01'
def run_model(dist, outcome_type="births", cat_name="total", rank=5, missingness=True, 
         disp_param=1e-4, sample_disp=False, placebo_state = None, 
         start_time = '2016-01-01', end_time = '2023-12-31',
         placebo_time = None, dobbs_donor_sensitivity=False,
         model_treated=True, results_file_suffix = "",
         num_chains=num_chains, num_warmup=1000, num_samples=1000, thinning=1):
    
    numpyro.set_host_device_count(num_chains)

    # df = pd.read_csv('data/dobbsbimonthlybirths_10_23_24.csv')
    df = pd.read_csv('/Users/shaokangyang/Library/CloudStorage/GoogleDrive-sky.ang510@gmail.com/My Drive/Code/dobbs_fertility/data/fertility_data_birth_update.csv')
    
    from clean_monthly_birth_data import prep_data, clean_dataframe, create_unit_placebo_dataset, create_time_placebo_dataset
    
    df = clean_dataframe(df, outcome_type, cat_name,  
                         dobbs_donor_sensitivity=dobbs_donor_sensitivity, csv_filename=None)
    df = df[df['time'] <= pd.to_datetime(end_time)]
    df = df.sort_values(by=['state', 'year', 'bmcode']) 
    df = df.drop_duplicates()

    if placebo_state is not None and placebo_state != "Texas":
        df = create_unit_placebo_dataset(df, placebo_state = placebo_state)
    
    if placebo_time is not None:
        df = create_time_placebo_dataset(df, new_treatment_start = placebo_time)
    else:
        # Only use data from 2016 onwards if not using a placebo time
        df = df[df['time'] >= pd.to_datetime(start_time)]  

    data_dict_cat = prep_data(df, outcome_type=outcome_type, group=cat_name)

    data_dict_cat['Y'].shape
    data_dict_cat['denominators'].shape
    data_dict_cat['control_idx_array'].shape
    
    import numpy as np
    from jax import random
    from numpyro.infer import MCMC, NUTS, Predictive

    #from models.monthly_model import monthly_model

    # set the random seed
    rng_key = random.PRNGKey(8675309)
    # split the random key
    rng_key, rng_key_ = random.split(rng_key)
    # Setup the sampler
    kernel = NUTS(model)

    mcmc = MCMC(
        kernel,
        num_warmup=num_warmup,
        num_samples=num_samples,
        num_chains=num_chains,
        progress_bar=True,
        thinning=thinning
    )

    mcmc.run(
        rng_key_,
        y=data_dict_cat['Y'],
        denominators=data_dict_cat['denominators'],
        control_idx_array=data_dict_cat['control_idx_array'],
        missing_idx_array=data_dict_cat['missing_idx_array'],
        rank=rank,
        outcome_dist=dist,
        adjust_for_missingness=missingness,
        nb_disp = disp_param,
        sample_disp = sample_disp,
        model_treated = model_treated
    )

    samples = mcmc.get_samples(group_by_chain=True)
    predictive = Predictive(model, mcmc.get_samples(group_by_chain=False))
    rng_key, rng_key_ = random.split(rng_key)

    predictions = predictive(
        rng_key_,
        denominators=data_dict_cat['denominators'],
        control_idx_array=None, #data_dict_cat['control_idx_array'],
        missing_idx_array=None, #data_dict_cat['missing_idx_array'],
        rank=rank,
        outcome_dist=dist,
        nb_disp = disp_param,
        sample_disp = sample_disp,
        model_treated = False
    )['y_obs']
    K, D, N = data_dict_cat['denominators'].shape
    pred_mat = predictions.reshape(mcmc.num_chains, int(mcmc.num_samples / mcmc.thinning), K, D, N)
   
    ## Take Python output and convert to draws matrix form
    params = dict_to_tidybayes({'mu': samples['mu_ctrl'], 'te': samples['te'], 'disp' : samples['disp']})
    preds = dict_to_tidybayes({"ypred" : pred_mat})

    preds[".chain"] = params[".chain"]
    preds[".draw"] = params[".draw"]

    all_samples = params.merge(preds, left_on = ['.draw', '.chain'], right_on = ['.draw', '.chain'])
    results_df = pd.DataFrame(all_samples)

    ## save input df
    df.to_csv('/Users/shaokangyang/Library/CloudStorage/GoogleDrive-sky.ang510@gmail.com/My Drive/Code/fertility_results/2024/df_{}.csv'.format(results_file_suffix))
    ## save posterior samples
    results_df.to_csv(
        '/Users/shaokangyang/Library/CloudStorage/GoogleDrive-sky.ang510@gmail.com/My Drive/Code/fertility_results/2024/{}_{}_{}_{}_{}.csv'.format(dist, "births", cat_name, rank, results_file_suffix)
    )

    
if __name__ == '__main__':
    from clean_monthly_birth_data import subgroup_definitions
    from joblib import Parallel, delayed

    # Define the inputs for the function
    inputs = [6, 7, 8, 9, 10, 11, 12]
    outcome_type = "births" 
    cats = list(subgroup_definitions[outcome_type].keys())
    dists = ['NB'] # Poisson or NB
    missing_flags = [True]
    # disp_params = [1e-4, 1e-3]
    disp_params = [1e-4]
    ## placebo_times = ["2020-05-01"]
    placebo_times = [None]
    placebo_states = [None]
    sample_disp = False
    dobbs_donor_sensitivity = False

    args = [(dist, cat, rank, m, disp, p, tm) for dist in dists for rank in inputs for cat in cats 
            for m in missing_flags for disp in disp_params for p in placebo_states 
            for tm in placebo_times]
    # Run the function in parallel
    results = Parallel(n_jobs=8)(delayed(run_model)(dist=i[0], outcome_type=outcome_type, cat_name=i[1], rank=i[2], missingness=i[3], 
                                                disp_param=i[4],
                                                sample_disp=sample_disp, placebo_state=i[5], placebo_time = i[6], 
                                                dobbs_donor_sensitivity=dobbs_donor_sensitivity, 
                                                results_file_suffix="through_june", num_chains=4, num_samples=2500, num_warmup=1000, thinning=10) for i in args)

