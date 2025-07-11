import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols
from sklearn.preprocessing import StandardScaler

def solve_task():
    """
    This function simulates the data generation and model comparison from the R script
    to determine the expected p-value from the anova test.
    """
    # Set a seed for reproducibility
    np.random.seed(42)
    n = 100

    # 1. Replicate the R data generation process
    # In R, `1:2` is recycled to the length of the dataframe (n=100)
    water_availability = np.tile([1, 2], n // 2)
    plant_intra_div = np.random.randint(1, 3, n)
    plant_inter_div = np.random.randint(1, 11, n)

    # Create the initial dataframe
    df = pd.DataFrame({
        'water_availability': water_availability,
        'plant_intra_div': plant_intra_div,
        'plant_inter_div': plant_inter_div
    })

    # Generate the dependent variables based on the specified formulas
    # The key is the strong effect (3.14) of insect_richness on herbivory
    df['insect_richness'] = (df['water_availability'] * 0.01 +
                             df['plant_intra_div'] * 0.5 +
                             df['plant_inter_div'] * 1.2 +
                             np.random.randn(n))

    df['herbivory'] = (df['insect_richness'] * 3.14 +
                       df['water_availability'] * 0.5 +
                       df['plant_intra_div'] * 0.1 +
                       df['plant_inter_div'] * 0.2 +
                       np.random.randn(n))

    # Replicate the `scale()` function from R
    scaler = StandardScaler()
    df_scaled = pd.DataFrame(scaler.fit_transform(df), columns=df.columns)

    # 2. Fit the two nested models
    # Model 1 (full model)
    model_1_formula = 'herbivory ~ insect_richness + water_availability + plant_intra_div + plant_inter_div'
    fit_1 = ols(model_1_formula, data=df_scaled).fit()

    # Model 2 (restricted model, without insect_richness)
    model_2_formula = 'herbivory ~ water_availability + plant_intra_div + plant_inter_div'
    fit_2 = ols(model_2_formula, data=df_scaled).fit()

    # 3. Perform ANOVA to compare the two nested models
    # This is the Python equivalent of R's `anova(fit_1, fit_2)` for this problem
    anova_results = sm.stats.anova_lm(fit_2, fit_1)

    # 4. Print the results
    # The anova_lm function returns a dataframe. We extract the key values.
    # The comparison is on the second row of the results table.
    df_diff = anova_results.loc[1, 'df_diff']
    f_statistic = anova_results.loc[1, 'F']
    p_value = anova_results.loc[1, 'Pr(>F)']

    print("ANOVA Comparison of Nested Models")
    print("="*35)
    print(f"The 'final equation' being tested is whether the fit of Model 1 is significantly better than Model 2.")
    print(f"This is tested with an F-test (equivalent to the Chi-square test in the R code).")
    print("\n--- Test Results ---")
    print(f"Difference in Degrees of Freedom: {df_diff:.0f}")
    print(f"F-statistic for the comparison: {f_statistic:.4f}")
    print(f"P-value (Pr(>F)): {p_value}")
    print("="*35)
    print("\nConclusion: The P-value is extremely small (effectively zero), indicating that the path from")
    print("'insect_richness' to 'herbivory' is highly significant, and Model 1 is a much better fit.")

solve_task()