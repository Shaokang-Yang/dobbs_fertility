import numpy as np
import pandas as pd
import re


births_subgroup_definitions = {
        'race': ("nhwhite", "hisp", "nhblack", "otherraceeth"),
        'edu': ("nohs", "hs", "somecoll", "coll"),
       # 'edu': ("hs_less", "somecoll_more"),
        #'age': ("age1519", "age2024", "age2529", "age3034", "age3539", "age4044"),
        'age': ("age1524", "age2534", "age3544"),
        'insurance': ("medicaid", "nonmedicaid"),
        # California should be dropped for marital. 
        'marital': ("married", "unmarried"),
        'total': ("total",),
    }
deaths_subgroup_definitions = {
        'race': ("nhwhite", "hisp", "nhblack", "otherraceeth"),
        'neonatal' : ("neo", "nonneo"),
        'congenital' : ("con", "noncon"),
        'total': ("total",),
    }

subgroup_definitions = {
    'births': births_subgroup_definitions,
    'deaths': deaths_subgroup_definitions
}

def clean_dataframe(dat:pd.DataFrame, outcome_type="births", cat_name="total", 
                    csv_filename='data/dat_quarterly.csv', drop_dobbs=False, 
                    dobbs_donor_sensitivity=False):
    """
    Filters, imputes, and adds relevant columns to the dataframe
    """

    # Set time variable for bimonths
    dat['time'] = pd.to_datetime(dat.year.astype(str) + '-' + (dat.bmcode * 2 - 1).astype(str) + "-01")

    dat['deaths_nonneo'] = dat['deaths_total'] - dat['deaths_neo']
    dat['deaths_noncon'] = dat['deaths_total'] - dat['deaths_con']
    dat['births_con'] = dat['births_noncon'] = dat['births_total']
    dat['births_neo'] = dat['births_nonneo'] = dat['births_total']
    
    # Convert "births_nhblack" column to numeric, replacing "Suppressed" with NaN
    dat['births_nhblack'] = pd.to_numeric(dat['births_nhblack'].replace("Suppressed", pd.NA))

    # Create a new column 'partial_ban' based on conditions from 'dobbscodev3' column
    dat['partial_ban'] = dat['dobbscodev3'].apply(lambda x: 1 if x == 2 else 0)
    #dat['time'] = pd.to_datetime(dat.year.astype(str) + '-' + dat.month, format="%Y-%B")

    def fill_in_missing_denoms(dat):
        # Get a list of column names containing the string "pop"
        cols_with_pop = dat.filter(regex=r'pop').columns
        print(dat['state'].unique())
        # Iterate over each column containing "pop"
        for col in cols_with_pop:
            # Find the row index of the maximum value for the current column in 2022
            pop_index_2022 = dat.loc[dat['year'] == 2022, col].idxmax()
            # Find the row index of the maximum value for the current column in 2021
            pop_index_2021 = dat.loc[dat['year'] == 2021, col].idxmax()
            
            # Impute missing values in the current column
            # For rows with missing values, replace with the squared value of the 2022 row divided by the 2021 row
            dat.loc[dat[col].isna(), col] = (dat.loc[pop_index_2022, col] ** 2) / dat.loc[pop_index_2021, col]
        
        # Return the modified DataFrame
        return dat

    # Hacky imputation
    # All of 2023 currently has a population of `NA`
    # We'll use a linear imputation of (pop 2022) * (pop 2022) / (pop 2021)
    # Well do so by:
    ## Group the original DataFrame by 'state'
    ## Apply the fill_in_missing_denoms function to each group
    ## Ungroup the result and reset the index
    dat = dat.groupby('state').apply(fill_in_missing_denoms).reset_index(drop=True)

    # create a column that is YYYY-QQ for indexing later
    #dat = dat.assign(time=dat['year'].astype(str) + '-' + ((dat['q'] - 1) * 3 + 1).astype(str))
    dat['quarter'] = dat['time'].apply(lambda x: f"{pd.Period(x, freq='Q').start_time.year}-Q{pd.Period(x, freq='Q').quarter}")


    ## Correct for different number of days each month

    # Assuming 'dat' is a DataFrame in Python
    # Convert 'time' column to datetime if it's not already
    dat['time'] = pd.to_datetime(dat['time'])

    # Filter the DataFrame to keep rows with time before 2024-12-01
    dat = dat[dat['time'] < pd.to_datetime('2024-01-01')]

    # Calculate the interval in days and divide by 91 to get the multiplier
    dat['days_multiplier'] = (dat['time'] + pd.DateOffset(months=3) - dat['time']).dt.days / 91


    ## States where abortion was banned in July 2022 include Alabama, Arkansas, Mississippi, Missouri, Oklahoma, South Dakota, 
    # Texas, West Virginia, and Wisconsin. Kentucky and Louisiana banned abortion in August and Idaho and Tennessee 
    # did so in September 2022. North Dakota banned abortion in April 2023. Ohio and South Carolina had six-week bans in place 
    # for the first two months after the Dobbs decision and Georgia had a six-week ban beginning in August 2022.  
    # States where abortion was banned in July 2022
    # july_banned_states = ['Alabama', 'Arkansas', 'Mississippi', 'Missouri', 'Oklahoma', 'South Dakota', 'Texas', 'West Virginia', 'Wisconsin']

    # # States where abortion was banned in August 2022
    # august_banned_states = ['Kentucky', 'Louisiana']

    # # States where abortion was banned in September 2022
    # september_banned_states = ['Idaho', 'Tennessee']

    # # States where abortion was banned in April 2023
    # april_banned_states = ['North Dakota']

    # # States with six-week bans in place for the first two months after the Dobbs decision
    # six_week_ban_states = ['Ohio', 'South Carolina']

    # # Georgia had a six-week ban beginning in August 2022
    # georgia_six_week_ban = ['Georgia']

    # # Combine all banned states
    # banned_states = (
    #     july_banned_states +
    #     august_banned_states +
    #     september_banned_states +
    #     april_banned_states +
    #     six_week_ban_states +
    #     georgia_six_week_ban
    # )


    # # Assign dobbscodev2 to dobbs_code
    # dat['dobbscodev2'] = dat.groupby('state')['dobbscodev2'].ffill()
    
    # dat['dobbs_code'] = dat['dobbscodev2']

    # # Set dobbs_code to 0 for Texas for times < 2022-04-01
    # dat.loc[(dat['state'] == 'Texas') & (dat['time'] < pd.to_datetime('2022-04-01')), 'dobbs_code'] = 0

    # # Set dobbs_code to 0 for all other states for times < 2023-01-01
    # dat.loc[(dat['state'] != 'Texas') & (dat['time'] < pd.to_datetime('2023-01-01')), 'dobbs_code'] = 0

    # # Set dobbs_code to 0 for Kentucky and Louisiana for times before 2023-02-01
    # dat.loc[(dat['state'].isin(['Kentucky', 'Louisiana'])) & (dat['time'] < pd.to_datetime('2023-02-01')), 'dobbs_code'] = 0
    
    # # Set dobbs_code to 0 for Idaho and Tennessee for times before 2023-03-01
    # dat.loc[(dat['state'].isin(['Idaho', 'Tennessee'])) & (dat['time'] < pd.to_datetime('2023-03-01')), 'dobbs_code'] = 0

    # # Set dobbs_code to 0 for Georgia
    # dat.loc[(dat['state'] == 'Georgia'), 'dobbs_code'] = 0
    # dat.loc[(dat['state'].isin(['Georgia'])) & (dat['time'] > pd.to_datetime('2023-04-01')), 'dobbs_code'] = 1
    


        # Create a control index array DataFrame
    if outcome_type == "births":
        dat['exposure_code'] = dat['exposed_births']
    if outcome_type == "deaths":
        dat['exposure_code'] = dat['exposed_infdeaths']
    # Convert to a list of unique states with dobbs_code == 1
    
    states_with_ban = dat.loc[dat['exposure_code'] == 1, 'state'].unique().tolist()
    print(states_with_ban)

    # Create a new column 'births_other' by subtracting births of non-Hispanic white, Hispanic, and non-Hispanic black from total births
    dat['births_other'] = dat['births_total'] - dat['births_nhwhite'] - dat['births_hisp'] - dat['births_nhblack']

    # Create a new column 'pop_other' by subtracting population of non-Hispanic white, Hispanic, and non-Hispanic black from total population
    dat['pop_other'] = dat['pop_total'] - dat['pop_nhwhite'] - dat['pop_hisp'] - dat['pop_nhblack']
    dat = dat.sort_values(['state', 'time'])

    # Remove California for marital
    if cat_name == "marital":
        dat = dat[dat["state"] != "California"]
    
    # Drop Vermont, Wyoming and Montana for mortality
    if outcome_type == "deaths":
        dat = dat[~dat["state"].isin(["Vermont", "Wyoming", "Montana"])]
    
    if drop_dobbs:
        dat = dat[~(dat["state"].isin(states_with_ban) & (dat["state"] != 'Texas'))]

    if dobbs_donor_sensitivity:
        sensitivity_states = dat[~dat["dobbscode_sensitivity"].isna()]['state'].unique()
        dat = dat[dat["state"].isin(sensitivity_states)]

    if csv_filename is not None:
        ## Save to csv so we don't have to do this every time 
        dat.to_csv(csv_filename)
    return dat


def prep_data(dat, group=None, outcome_type="births", variables=None, covariates=None):
    """
    Prepare data for analysis by creating DataFrames for births or deaths (numerators), population or births (denominators), control indices, and missing indices.

    Args:
        dat (pandas.DataFrame): Input data containing information about births, population, and other relevant variables.
        variables (list, optional): List of variable names (e.g., "white", "hisp", "black", "other"). Default is ["white", "hisp", "black", "other"].
        covariates (list, optional): List of covariate names to include in the analysis. Default is None.

    Returns:
        dict: A dictionary containing the following items:
            Y (pandas.DataFrame): DataFrame with birth counts for each category, state, and time period.
            population (pandas.DataFrame): DataFrame with population counts for each category, state, and time period.
            state_fe (numpy.ndarray): Array of state fixed effects.
            control_idx_array (pandas.DataFrame): DataFrame with control indices for each category, state, and time period.
            missing_idx_array (pandas.DataFrame): DataFrame with missing indices for each category, state, and time period.
            days_multiplier (float): Days multiplier value.
            variables (list): List of variable names used in the analysis.
            D_cov (numpy.ndarray or None): Matrix of covariates, if provided. If no covariates are provided, D_cov is set to None.
    """
    if (group is not None) and (variables is not None):
        raise Exception("Only one of group/variables can be specified.")
    if group is not None:
        variables = subgroup_definitions[outcome_type][group]
    

    if outcome_type == "deaths":
        
        # Create a list of death column names
        death_columns = ["deaths_" + var for var in variables]    
        
        # Create a deaths DataFrame
        deaths = (
            dat[["state", "time"] + death_columns]  # Select 'state', 'time', and death columns
            .melt(id_vars=["state", "time"], value_vars=death_columns, var_name="category", value_name="deaths")  # Melt death columns into long format
            .pivot_table(index=["category"], columns=["state", "time"], values="deaths", aggfunc="sum", fill_value=0)  # Pivot to wide format, summing death values
        )

    # Create a list of birth column names
    birth_columns = ["births_" + var for var in variables]
    
    # Create a list of population column names
    denom_columns = ["pop_" + var for var in variables]

    if outcome_type == "births":
        # Create a population DataFrame
        population = (
            dat[["state", "time"] + denom_columns]  # Select 'state', 'time', and population columns
            .melt(id_vars=["state", "time"], value_vars=denom_columns, var_name="category", value_name="population")  # Melt population columns into long format
            .pivot_table(index=["category"], columns=["state", "time"], values="population", aggfunc="sum", fill_value=0)  # Pivot to wide format, summing population values 
        ) 

    # Create a births DataFrame
    births = (
        dat[["state", "time"] + birth_columns]  # Select 'state', 'time', and birth columns
        .melt(id_vars=["state", "time"], value_vars=birth_columns, var_name="category", value_name="births")  # Melt birth columns into long format
        .pivot_table(index=["category"], columns=["state", "time"], values="births", aggfunc="sum", fill_value=0)  # Pivot to wide format, summing birth values
    )
    
    if outcome_type == "deaths":
        Y = deaths
        denominators = births
        outcome_columns = death_columns
    else:
        Y = births
        denominators = population / 1e4 # Population per 1000
        outcome_columns = birth_columns
    
    control_idx_array = (
        dat[["state", "time", "exposure_code"] + outcome_columns]  # Select 'state', 'time', 'exposed_births', and birth/death columns
        .melt(id_vars=["state", "time", "exposure_code"], value_vars=outcome_columns, var_name="category", value_name=outcome_type)  # Melt birth columns into long format
        .assign(ctrl_index=(lambda x: x["exposure_code"] == 0))  # Create a control index column based on 'exposed_births'
        .pivot_table(index=["category"], columns=["state", "time"], values="ctrl_index", aggfunc="sum", fill_value=0)  # Pivot to wide format, summing control index values
    ).astype(np.bool_) # cast to a boolean so we don't have issues when we mask
    
    # Create a missing index array DataFrame
    missing_idx_array = (
        dat[["state", "time", "exposure_code"] + outcome_columns]  # Select 'state', 'time', 'exposure_code', and birth columns
        .melt(id_vars=["state", "time", "exposure_code"], value_vars=outcome_columns, var_name="category", value_name=outcome_type)  # Melt birth columns into long format
        .assign(missing_index=lambda x: x[outcome_type].isna().astype(int))  # Create a missing index column based on missing birth values
        .pivot_table(index=["category"], columns=["state", "time"], values="missing_index", aggfunc="sum", fill_value=0)  # Pivot to wide format, summing missing index values
    ).astype(np.bool_) # cast to a boolean so we don't have issues when we mask
    
    num_states = len(dat.state.unique())
    total_length = denominators.shape[1]
    denominators = denominators.values.reshape((len(variables), num_states, denominators.shape[1]//num_states))
    Y = Y.values.reshape((len(variables), num_states, total_length//num_states))

    # If covariates are provided, calculate the covariates matrix
    if covariates is not None:
        D_cov = dat.groupby("state")[covariates].mean().reset_index()[covariates].values
        D_cov[np.isnan(D_cov)] = D_cov[~np.isnan(D_cov)].mean()
    else:
        D_cov = None
    

    # Return a dictionary with the calculated values
    return {
        "Y": Y,
        "denominators": denominators,
        #"state_fe": state_fe,
        "control_idx_array": control_idx_array.values.reshape((len(variables), num_states, total_length//num_states)),
        "missing_idx_array": missing_idx_array.values.reshape((len(variables), num_states, total_length//num_states)),
        #"days_multiplier": days_multiplier,
        "variables": variables,
        "D_cov": D_cov,
    }

def create_unit_placebo_dataset(df, treated_state = "Texas", placebo_state = "California"):
    """
    Create a placebo dataset for by giving Texas' treatment times to `placebo_state` state by removing the state 
    and removing Texas.

    Args:
        df (pandas.DataFrame): Input data containing information about births, population, exposure codes and other relevant variables.
        placebo_state (str): Name of the placebo state to remove from the dataset.

    Returns:
        pandas.DataFrame: Placebo dataset with the specified placebo state with treatment times of treated_state.
    """
    
    print("Creating unit-placebo dataset for {}".format(placebo_state))

    # Get the columns that start with 'exposed' (exposed_births/exposed_deaths)
    exposure_columns = df.filter(regex='^exposure').columns

    # Get the values from the rows where 'state' equals 'treated_state'
    treated_values = df.loc[df['state'] == treated_state, exposure_columns]

    # Set the values in the rows where 'state' equals 'placebo_state' to the values from the 'treated_state' rows
    df.loc[df['state'] == placebo_state, exposure_columns] = treated_values.values
    
    # Filter the DataFrame to keep rows where 'state' is not equal to the treated state
    return df[df["state"] != treated_state]

def create_time_placebo_dataset(df, new_treatment_start="2022-05-01"):
    """
    Create a time placebo dataset by shifting treatment times early and capping the end date.
    
    Args:
        df (pandas.DataFrame): Input data containing information about births, population, exposure codes and other relevant variables.
        first_treatment_start (str): Start time of the first treated unit (usually Texas). Set to "2022-05-01" by default which is the actual SB8 time.

    Returns:
        pandas.DataFrame: Modified dataset with shifted time variables.
    """
    
    print("Creating placebo-in-time dataset starting in {}".format(new_treatment_start))
    
    # Convert 'first_treatment_start' to datetime
        
    new_treatment_start = pd.to_datetime(new_treatment_start)
    original_treatment_start = df.loc[df['exposure_code'] == 1, 'time'].min()

    end_date = df.loc[df['exposure_code'] == 1, 'time'].max()
    
    new_end = new_treatment_start + (end_date - original_treatment_start)
    
    original_time_length = end_date - pd.to_datetime("2016-01-01")
    new_start = new_end - original_time_length
    if new_start < df["time"].min():
        new_start = df["time"].min()
    
    new_start = new_start.to_period('M').to_timestamp()
    new_end = new_end.to_period('M').to_timestamp()
    
    new_time_length = new_end - new_start

    # Get the columns that start with 'exposure_code'
    exposure_code_values = df.loc[(df['time'] >= end_date - new_time_length), ['exposure_code']]

    df = df[(df['time'] <= new_end) & (df['time'] >= new_start)]

    print(exposure_code_values.values.shape)
    print(df.shape)
    if len(exposure_code_values) == len(df):
        df.loc[:, "exposure_code"] = exposure_code_values.values
    else:
        raise ValueError("The length of new exposure_code values does not match the number of rows in df")

    return df

if __name__ == '__main__':
    import argparse
    import pandas as pd
    import pickle 
    import os
    import gzip 
    
    parser = argparse.ArgumentParser(description='Script for cleaning and parsing quarterly births. We will assume that the file name is in the format quarterly_fertility_mortality_MMDDYY.csv where MM is month,DD is day, and YY is year.')
    parser.add_argument('filename', type=str, help='Name of the input file')
    parser.add_argument('--save-dict', action='store_true', help='Create and save the dictionary used for model training as a pkl file')
    parser.add_argument("--group", help="Subgroups to create a dictionary for", default="total")

    args = parser.parse_args()
    # Create a directory to hold the parsed files
    directory_name = 'model_data/' + args.filename.split('/')[-1].split('.')[0]
    os.makedirs(directory_name, exist_ok=True)
    dat = clean_dataframe(pd.read_csv(args.filename))
    # Save the dataframe
    dat.to_csv(directory_name + '/data_frame.csv')

    if args.save_dict:
        data_dict = prep_data(dat, subgroup_definitions[outcome_type][args.group])
        # Specify the filename for the gzipped pickle file
        dict_filename = 'births_{}_dict.pkl.gz'.format(args.group)

        # Open the file in binary write mode using gzip
        with gzip.open('{}/{}'.format(directory_name, dict_filename), 'wb') as f:
            # Use pickle to dump the data dictionary into the file
            pickle.dump(data_dict, f)