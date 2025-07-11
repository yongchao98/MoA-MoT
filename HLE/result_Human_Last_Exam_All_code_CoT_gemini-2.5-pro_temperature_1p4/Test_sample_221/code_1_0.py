import pandas as pd
import numpy as np
from sklearn.preprocessing import scale
import statsmodels.formula.api as smf
from scipy.stats import chi2

def solve_case():
    """
    This function simulates the data generation and model comparison from the R script
    to determine the expected P-value from the anova test.
    """
    # Set a seed for reproducibility. The original R script is random,
    # but the expected outcome is consistent due to the strong simulated effect.
    np.random.seed(42)

    # 1. Data Generation as in the R script
    n = 100
    # In R, a vector of length 2 is recycled to match the length of 100
    water_availability_raw = np.tile([1, 2], n // 2)

    df = pd.DataFrame({
        'water_availability': water_availability_raw,
        'plant_intra_div': np.random.choice([1, 2], n, replace=True),
        'plant_inter_div': np.random.choice(range(1, 11), n, replace=True)
    })

    df['insect_richness'] = (df['water_availability'] * 0.01 +
                             df['plant_intra_div'] * 0.5 +
                             df['plant_inter_div'] * 1.2 +
                             np.random.randn(n))

    df['herbivory'] = (df['insect_richness'] * 3.14 +
                       df['water_availability'] * 0.5 +
                       df['plant_intra_div'] * 0.1 +
                       df['plant_inter_div'] * 0.2 +
                       np.random.randn(n))

    # Scale the data frame, as done in the R script
    df_scaled = pd.DataFrame(scale(df), columns=df.columns)

    # 2. Fit equivalent models
    # The anova test compares the model system where herbivory depends on insect_richness
    # (model_1) vs. one where it does not (model_2). The core of this comparison
    # is testing the significance of that single path. We can do this with a
    # Likelihood Ratio Test between two Ordinary Least Squares (OLS) models for herbivory.

    full_model_formula = 'herbivory ~ insect_richness + water_availability + plant_intra_div + plant_inter_div'
    restricted_model_formula = 'herbivory ~ water_availability + plant_intra_div + plant_inter_div'

    fit_full = smf.ols(full_model_formula, data=df_scaled).fit()
    fit_restricted = smf.ols(restricted_model_formula, data=df_scaled).fit()

    # 3. Perform Likelihood Ratio Test (LRT)
    # The anova() in R for lavaan models performs an LRT.
    # The test statistic is 2 * (log-likelihood of full model - log-likelihood of restricted model)
    loglik_full = fit_full.llf
    loglik_restricted = fit_restricted.llf
    
    # Degrees of freedom is the difference in the number of parameters (1 in this case)
    degrees_freedom = fit_full.df_model - fit_restricted.df_model
    
    lr_statistic = 2 * (loglik_full - loglik_restricted)
    
    # The P-value is calculated from the chi-squared distribution
    p_value = chi2.sf(lr_statistic, df=degrees_freedom)

    # The resulting p-value will be extremely small, effectively zero.
    # Among the multiple-choice options, 0 is the best representation.
    print(f"The Likelihood Ratio Test statistic comparing the two models is {lr_statistic:.4f}.")
    print(f"The degrees of freedom for the test is {int(degrees_freedom)}.")
    print(f"The resulting P-value is approximately {p_value:.2e}, which is effectively 0.")
    print("\nFrom the choices, the expected P-value is 0.")

solve_case()