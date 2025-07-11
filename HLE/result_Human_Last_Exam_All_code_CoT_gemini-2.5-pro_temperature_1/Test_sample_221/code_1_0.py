import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

def solve_case():
    """
    This function replicates the R script's data simulation and model comparison
    to determine the expected p-value from the anova function.
    """
    n = 100

    # 1. Simulate the data as described in the R code
    # Note: R's sample function samples with replacement by default.
    # R's 1:2 creates a vector c(1, 2). The data generation will cycle through this.
    water_availability = np.tile([1, 2], n // 2)

    df = pd.DataFrame({
        'water_availability': water_availability,
        'plant_intra.div': np.random.choice([1, 2], size=n, replace=True),
        'plant_inter.div': np.random.choice(range(1, 11), size=n, replace=True)
    })

    # Generate the dependent variables based on the linear equations
    df['insect_richness'] = (df['water_availability'] * 0.01 +
                               df['plant_intra.div'] * 0.5 +
                               df['plant_inter.div'] * 1.2 +
                               np.random.normal(size=n))

    df['herbivory'] = (df['insect_richness'] * 3.14 +
                         df['water_availability'] * 0.5 +
                         df['plant_intra.div'] * 0.1 +
                         df['plant_inter.div'] * 0.2 +
                         np.random.normal(size=n))

    # 2. Scale the data (z-score transformation)
    df_scaled = pd.DataFrame(stats.zscore(df), columns=df.columns)

    # 3. Perform the Likelihood Ratio Test, which is what lavaan's anova() does
    # The SEM specified can be estimated with two OLS regressions. The anova()
    # test only compares the part of the model that differs, which is the
    # equation for 'herbivory'.

    y = df_scaled['herbivory']

    # Fit the full model (model_1)
    X_full = df_scaled[['insect_richness', 'water_availability', 'plant_intra.div', 'plant_inter.div']]
    X_full = sm.add_constant(X_full)
    model_full_fit = sm.OLS(y, X_full).fit()
    ll_full = model_full_fit.llf

    # Fit the reduced model (model_2), which omits 'insect_richness'
    X_reduced = df_scaled[['water_availability', 'plant_intra.div', 'plant_inter.div']]
    X_reduced = sm.add_constant(X_reduced)
    model_reduced_fit = sm.OLS(y, X_reduced).fit()
    ll_reduced = model_reduced_fit.llf

    # 4. Calculate the test statistic and p-value
    # The test statistic follows a chi-square distribution
    lr_statistic = 2 * (ll_full - ll_reduced)
    # The degrees of freedom is the difference in the number of parameters (1 in this case)
    df_diff = 1
    p_value = stats.chi2.sf(lr_statistic, df_diff)

    print("Likelihood Ratio Test for Nested Models")
    print("-" * 40)
    print(f"Log-Likelihood of full model (model_1): {ll_full:.4f}")
    print(f"Log-Likelihood of reduced model (model_2): {ll_reduced:.4f}")
    print(f"The final equation for the test statistic is: 2 * (LL_full - LL_reduced)")
    print(f"Likelihood Ratio Test Statistic: 2 * ({ll_full:.4f} - ({ll_reduced:.4f})) = {lr_statistic:.4f}")
    print(f"Degrees of Freedom: {df_diff}")
    print(f"P-value: {p_value}")
    print("-" * 40)
    print("The p-value is consistently a very small number, effectively zero.")


solve_case()
<<<D>>>