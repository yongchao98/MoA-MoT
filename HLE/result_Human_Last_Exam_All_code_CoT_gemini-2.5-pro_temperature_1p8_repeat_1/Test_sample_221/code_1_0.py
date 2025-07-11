import numpy as np
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
import pandas as pd

def solve():
    """
    Simulates the data and model comparison described in the R script to find the p-value.
    """
    # Set a seed for reproducibility, similar to R's set.seed()
    np.random.seed(42)

    n = 100

    # 1. Generate data
    # In R, water_availability = 1:2. In a dataframe of n=100, this becomes a repeating sequence.
    water_availability = np.tile([1, 2], n // 2)

    df = pd.DataFrame({
        'water_availability': water_availability,
        'plant_intra_div': np.random.randint(1, 3, n),
        'plant_inter_div': np.random.randint(1, 11, n)
    })

    # Generate dependent variables based on the formulas
    df['insect_richness'] = (df['water_availability'] * 0.01 +
                             df['plant_intra_div'] * 0.5 +
                             df['plant_inter_div'] * 1.2 +
                             np.random.randn(n))

    df['herbivory'] = (df['insect_richness'] * 3.14 +
                       df['water_availability'] * 0.5 +
                       df['plant_intra_div'] * 0.1 +
                       df['plant_inter_div'] * 0.2 +
                       np.random.randn(n))

    # 2. Scale the data (z-score normalization)
    df_scaled = pd.DataFrame(stats.zscore(df), columns=df.columns)

    # 3. Fit the two models
    # The anova comparison in lavaan tests the difference between the herbivory equations.
    
    # Model 1: Includes insect_richness
    model_1_formula = 'herbivory ~ insect_richness + water_availability + plant_intra_div + plant_inter_div'
    fit_1 = smf.ols(model_1_formula, data=df_scaled).fit()

    # Model 2: Excludes insect_richness
    model_2_formula = 'herbivory ~ water_availability + plant_intra_div + plant_inter_div'
    fit_2 = smf.ols(model_2_formula, data=df_scaled).fit()

    # 4. Perform the Likelihood Ratio Test (equivalent to lavaan's anova)
    lr_statistic = 2 * (fit_1.llf - fit_2.llf)
    # The degrees of freedom is the difference in the number of parameters (1 in this case)
    df_diff = fit_1.df_model - fit_2.df_model
    
    # The p-value is calculated from the chi-squared distribution
    p_value = stats.chi2.sf(lr_statistic, df=df_diff)

    # Print the equation representing the test
    # Note: Values are rounded for presentation.
    # Chisq = 2 * (logLik(model1) - logLik(model2))
    print(f"Likelihood-Ratio Test:")
    print(f"Chi-squared statistic = 2 * (LogLik_model_1 - LogLik_model_2)")
    print(f"{lr_statistic:.4f} = 2 * ({fit_1.llf:.4f} - ({fit_2.llf:.4f}))")
    print(f"Degrees of freedom = {int(df_diff)}")
    print(f"P-value = {p_value}")


solve()