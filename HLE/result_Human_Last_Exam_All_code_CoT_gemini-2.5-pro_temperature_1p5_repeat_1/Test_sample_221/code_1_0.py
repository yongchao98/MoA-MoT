import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm

def solve():
    """
    This function simulates the data generating process from the R script and
    performs a likelihood ratio test to determine the p-value from comparing
    the two models, which is what lavaan's anova() function does for nested models.
    """
    # Set a seed for reproducibility
    np.random.seed(42)
    n = 100

    # 1. Simulate the data as in the R script
    # R's recycling rule is emulated with np.tile
    example_df = pd.DataFrame({
        'water_availability': np.tile([1, 2], n // 2),
        'plant_intra.div': np.random.choice([1, 2], size=n, replace=True),
        'plant_inter.div': np.random.choice(np.arange(1, 11), size=n, replace=True)
    })

    # Generate the dependent variables based on the formulas
    # rnorm(n) is equivalent to np.random.normal(size=n)
    example_df['insect_richness'] = (example_df['water_availability'] * 0.01 +
                                  example_df['plant_intra.div'] * 0.5 +
                                  example_df['plant_inter.div'] * 1.2 +
                                  np.random.normal(size=n))

    example_df['herbivory'] = (example_df['insect_richness'] * 3.14 +
                            example_df['water_availability'] * 0.5 +
                            example_df['plant_intra.div'] * 0.1 +
                            example_df['plant_inter.div'] * 0.2 +
                            np.random.normal(size=n))

    # Scale the data, same as R's scale()
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(example_df)
    example_scaled = pd.DataFrame(scaled_data, columns=example_df.columns)

    # 2. Define and fit the two models for the 'herbivory' equation
    # The anova() test focuses on the difference between the herbivory equations.
    # This is equivalent to a likelihood ratio test between two linear models.

    Y = example_scaled['herbivory']

    # Model 1 (full model) includes insect_richness
    X_full = example_scaled[['insect_richness', 'water_availability', 'plant_intra.div', 'plant_inter.div']]
    X_full = sm.add_constant(X_full)  # Add intercept
    model_full_fit = sm.OLS(Y, X_full).fit()
    ll_full = model_full_fit.llf

    # Model 2 (reduced model) omits insect_richness
    X_reduced = example_scaled[['water_availability', 'plant_intra.div', 'plant_inter.div']]
    X_reduced = sm.add_constant(X_reduced)  # Add intercept
    model_reduced_fit = sm.OLS(Y, X_reduced).fit()
    ll_reduced = model_reduced_fit.llf

    # 3. Perform the Likelihood Ratio (chi-square difference) test
    # The test statistic
    lr_statistic = 2 * (ll_full - ll_reduced)

    # Degrees of freedom is the difference in the number of parameters
    df_difference = model_full_fit.df_model - model_reduced_fit.df_model

    # The p-value from the chi-squared distribution
    p_value = stats.chi2.sf(lr_statistic, df=df_difference)
    
    # 4. Print the results showing all the numbers in the calculation
    print("Likelihood Ratio Test (equivalent to lavaan's anova comparison)")
    print(f"Chi-square statistic calculation:")
    print(f"Statistic = 2 * (LogLikelihood_Model1 - LogLikelihood_Model2)")
    print(f"Statistic = 2 * ({ll_full:.2f} - ({ll_reduced:.2f}))")
    print(f"Calculated Statistic = {lr_statistic:.4f}")
    print(f"\nDegrees of Freedom = (Params in Model 1) - (Params in Model 2)")
    print(f"Degrees of Freedom = {int(model_full_fit.df_model)} - {int(model_reduced_fit.df_model)} = {int(df_difference)}")
    print(f"\nFinal P-value from Chi-square test: {p_value}")
    print(f"\nThis P-value is effectively zero.")


solve()