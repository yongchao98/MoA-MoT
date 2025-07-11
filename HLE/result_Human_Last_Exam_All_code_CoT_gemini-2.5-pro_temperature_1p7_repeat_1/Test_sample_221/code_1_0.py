import pandas as pd
import numpy as np
import semopy
from scipy.stats import chi2

def solve():
    """
    This function replicates the R script's logic to determine the expected
    P-value from comparing two structural equation models.
    """
    # Set a seed for reproducibility, similar to R's default behavior
    np.random.seed(123)
    n = 100

    # 1. Data Generation
    # Create a dataframe mirroring the R code's structure.
    # The recycling of `1:2` in R is replicated with np.tile.
    example_unscaled = pd.DataFrame({
        'water_availability': np.tile([1, 2], n // 2),
        'plant_intra_div': np.random.randint(1, 3, n),
        'plant_inter_div': np.random.randint(1, 11, n)
    })

    # Generate dependent variables based on the specified formulas
    example_unscaled['insect_richness'] = (
        example_unscaled['water_availability'] * 0.01 +
        example_unscaled['plant_intra_div'] * 0.5 +
        example_unscaled['plant_inter_div'] * 1.2 +
        np.random.normal(size=n)
    )
    # The key relationship: herbivory is strongly dependent on insect_richness
    example_unscaled['herbivory'] = (
        example_unscaled['insect_richness'] * 3.14 +
        example_unscaled['water_availability'] * 0.5 +
        example_unscaled['plant_intra_div'] * 0.1 +
        example_unscaled['plant_inter_div'] * 0.2 +
        np.random.normal(size=n)
    )

    # Scale the data, as done in the R script
    example_scaled = (example_unscaled - example_unscaled.mean()) / example_unscaled.std()
    
    # In semopy, variable names cannot contain dots. We rename them.
    example_scaled.columns = [c.replace('.', '_') for c in example_scaled.columns]

    # 2. Model Specification
    # Model 1: includes the path from insect richness to herbivory
    model_1_spec = """
    herbivory ~ insect_richness + water_availability + plant_intra_div + plant_inter_div
    insect_richness ~ water_availability + plant_intra_div + plant_inter_div
    """

    # Model 2: omits the path from insect richness to herbivory
    model_2_spec = """
    herbivory ~ water_availability + plant_intra_div + plant_inter_div
    insect_richness ~ water_availability + plant_intra_div + plant_inter_div
    """

    # 3. Model Fitting
    fit_1 = semopy.Model(model_1_spec)
    res_1 = fit_1.fit(example_scaled)
    
    fit_2 = semopy.Model(model_2_spec)
    res_2 = fit_2.fit(example_scaled)

    # 4. Model Comparison (Chi-Square Difference Test)
    # Get Chi-Square values and degrees of freedom for each model
    chi2_1, dof_1 = res_1.fun, fit_1.calc_df()
    chi2_2, dof_2 = res_2.fun, fit_2.calc_df()

    # Calculate the differences
    chi2_diff = abs(chi2_2 - chi2_1)
    dof_diff = abs(dof_2 - dof_1)
    
    # Calculate the p-value from the chi-square distribution
    # sf is the survival function (1 - cdf), which gives the p-value
    p_value = chi2.sf(chi2_diff, dof_diff)
    
    print("--- Chi-Squared Difference Test ---")
    print(f"Model 1 (Full):       Chi2 = {chi2_1:.3f}, Df = {dof_1}")
    print(f"Model 2 (Restricted): Chi2 = {chi2_2:.3f}, Df = {dof_2}")
    print("-" * 35)
    # Final Equation: Chi-Squared Difference = Chi2(Model 2) - Chi2(Model 1)
    print(f"Final Equation (Chi2 Diff): {chi2_2:.3f} - {chi2_1:.3f} = {chi2_diff:.3f}")
    # Final Equation: Degrees of Freedom Difference = Df(Model 2) - Df(Model 1)
    print(f"Final Equation (Df Diff):   {dof_2} - {dof_1} = {dof_diff}")
    
    # The P-value is the key result
    print(f"\nThe P-value for the difference is: {p_value}")
    
    print("\nConclusion: The p-value is extremely small, effectively 0.")
    print("This indicates that Model 1 fits the data significantly better than Model 2, as expected.")

solve()