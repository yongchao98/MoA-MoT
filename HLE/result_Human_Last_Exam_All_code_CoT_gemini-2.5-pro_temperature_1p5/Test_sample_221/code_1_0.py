import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.formula.api as smf

def solve():
    """
    Simulates the dataset and performs a likelihood ratio test to find
    the p-value for comparing the two nested models.
    """
    # Set a seed for reproducibility, so the random data is the same each time
    np.random.seed(123)
    n = 100

    # Generate data in a way that mimics the R script
    # R's 1:2 is recycled to fill n=100 observations: 1, 2, 1, 2, ...
    water_availability = np.tile([1, 2], n // 2)

    # Create the initial dataframe
    df_raw = pd.DataFrame({
        'water_availability': water_availability,
        'plant_intra_div': np.random.randint(1, 3, n),
        'plant_inter_div': np.random.randint(1, 11, n)
    })

    # Create the dependent variables based on the specified formulas
    df_raw['insect_richness'] = (df_raw['water_availability'] * 0.01 +
                                 df_raw['plant_intra_div'] * 0.5 +
                                 df_raw['plant_inter_div'] * 1.2 +
                                 np.random.randn(n))

    df_raw['herbivory'] = (df_raw['insect_richness'] * 3.14 +
                           df_raw['water_availability'] * 0.5 +
                           df_raw['plant_intra_div'] * 0.1 +
                           df_raw['plant_inter_div'] * 0.2 +
                           np.random.randn(n))

    # Scale the data using z-score, equivalent to R's scale()
    df_scaled = pd.DataFrame(stats.zscore(df_raw), columns=df_raw.columns)

    # The R anova() performs a likelihood ratio test (chi-square difference test).
    # We can replicate this by fitting two models for 'herbivory' and comparing them.

    # Model 1 (full model) for the herbivory equation, matches model_1
    fit_1_herbivory = smf.ols(
        'herbivory ~ insect_richness + water_availability + plant_intra_div + plant_inter_div',
        data=df_scaled
    ).fit()

    # Model 2 (reduced model) for the herbivory equation, matches model_2
    fit_2_herbivory = smf.ols(
        'herbivory ~ water_availability + plant_intra_div + plant_inter_div',
        data=df_scaled
    ).fit()

    # The difference in degrees of freedom is 1 (the single omitted parameter)
    df_diff = fit_1_herbivory.df_model - fit_2_herbivory.df_model

    # The test statistic is the Chi-squared difference from log-likelihoods
    chi_sq_diff = 2 * (fit_1_herbivory.llf - fit_2_herbivory.llf)

    # The p-value is calculated from the chi-squared distribution
    p_value = stats.chi2.sf(chi_sq_diff, df_diff)

    # Print the "equation" results from the anova test
    print("Likelihood Ratio Test (Chi-squared difference test) Results:")
    print(f"Chi-squared difference = {chi_sq_diff:.4f}")
    print(f"Degrees of freedom difference = {int(df_diff)}")
    # The p-value will be extremely small, printing as 0.0
    print(f"P-value (Pr(>Chisq)) = {p_value}")

solve()