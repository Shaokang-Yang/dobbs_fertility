import numpy as np
import pandas as pd

def dict_to_tidybayes(samples_dict, output_file=None):
    """
    Convert a dictionary of MCMC samples to a tidybayes draws dataframe.

    Parameters:
    - samples_dict: Dictionary with keys as parameter names and values as MCMC samples

    Returns:
    - draws_df: Tidybayes draws dataframe
    """
    # Get parameter names and dimensions
    param_names = list(samples_dict.keys())
    print(param_names)
    param_dims = [samples_dict[param].shape[2:] for param in param_names]

    vals = {}
    # Extract chains, draws
    chains, draws = samples_dict[param_names[0]].shape[:2]
    # Add columns for the chain and draw number
    vals['.chain'] = np.tile(np.arange(chains), draws)
    vals['.draw'] = np.repeat(np.arange(draws), chains)
    # Extract the samples 
    for param_name, param_dim in zip(param_names, param_dims):
        param_samples = samples_dict[param_name]
    
        # Reshape the samples to a 2D array (draws * chains, parameters)
        samples_array = param_samples.reshape((draws * chains, int(np.prod(param_dim))), order='F')
        # Create columns for each dimension of the parameter
        if len(param_dim) > 0:
            # Generate indices as strings
            index_strs = (np.array(np.meshgrid(*[np.arange(dim) for dim in param_dim], indexing='ij')).reshape(len(param_dim), np.prod(param_dim), order='F').T)
            for idx, index in enumerate(index_strs):
                vals[f"{param_name}[{','.join(map(str, index.flatten()+1))}]"] = samples_array[:, idx]
            
        else:
            vals[param_name] = samples_array[:, 0]
    
    draws_df = pd.DataFrame(vals)

    if output_file is not None:
        draws_df.to_csv(output_file, index_label='index')

    return draws_df