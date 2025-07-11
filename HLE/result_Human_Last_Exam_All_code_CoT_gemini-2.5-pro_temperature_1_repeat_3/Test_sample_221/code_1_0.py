import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from semopy import Model
from scipy.stats import chi2

def solve_task():
    """
    This function simulates the dataset and performs a likelihood ratio test
    to find the p-value from comparing the two nested SEM models, as described
    in the R script.
    """
    # Set a seed for reproducibility
    np.random.seed(42)

    # 1. Replicate Data Generation
    n = 100
    # In R, 1:2 is recycled to length n. np.tile accomplishes this.
    water_availability = np.tile([1, 2], n // 2)
    # In R, sample(1:2, n, replace=T) is equivalent to randint(1, 3, size=n)
    plant_intra_div = np.random.randint(1, 3, size=n)
    # In R, sample(1:10, n, replace=T) is equivalent to randint(1, 11, size=n)
    plant_inter_div = np.random.randint(1, 11, size=n)

    example_unscaled = pd.DataFrame({
        'water_availability': water_availability,
        'plant_intra_div': plant_intra_div,
        'plant_inter_div': plant_inter_div
    })

    # Generate dependent variables based on the R code logic
    # Using .loc to avoid potential SettingWithCopyWarning
    example_unscaled.loc[:, 'insect_richness'] = (
        example_unscaled['water_availability'] * 0.01 +
        example_unscaled['plant_intra_div'] * 0.5 +
        example_unscaled['plant_inter_div'] * 1.2 +
        np.random.randn(n)
    )
    example_unscaled.loc[:, 'herbivory'] = (
        example_unscaled['insect_richness'] * 3.14 +
        example_unscaled['water_availability'] * 0.5 +
        example_unscaled['plant_intra_div'] * 0.1 +
        example_unscaled['plant_inter_div'] * 0.2 +
        np.random.randn(n)
    )

    # Scale the data as in the R script
    scaler = StandardScaler()
    example_scaled = scaler.fit_transform(example_unscaled)
    # The original R code uses '.' in names, which semopy doesn't like. We replace with '_'.
    example = pd.DataFrame(example_scaled, columns=[c.replace('.', '_') for c in example_unscaled.columns])

    # 2. Define and Fit SEM Models
    model_1_spec = """
    herbivory ~ insect_richness + water_availability + plant_intra_div + plant_inter_div
    insect_richness ~ water_availability + plant_intra_div + plant_inter_div
    """

    model_2_spec = """
    herbivory ~ water_availability + plant_intra_div + plant_inter_div
    insect_richness ~ water_availability + plant_intra_div + plant_inter_div
    """

    fit_1 = Model(model_1_spec)
    fit_1.fit(example)

    fit_2 = Model(model_2_spec)
    fit_2.fit(example)

    # 3. Perform Likelihood Ratio Test (Chi-Square Difference Test)
    # The test statistic is the difference in the model chi-squares (-2*logLikelihood).
    # semopy stores this value in the 'objective_' attribute after fitting.
    chi_square_model_1 = fit_1.objective_
    chi_square_model_2 = fit_2.objective_
    
    test_statistic = chi_square_model_2 - chi_square_model_1

    # The degrees of freedom is the difference in the number of free parameters.
    # Model 1 has one more parameter (the path from insect_richness to herbivory) than Model 2.
    df_diff = 1
    
    # Calculate the p-value from the chi-square distribution
    p_value = chi2.sf(test_statistic, df_diff)

    # 4. Print the results, showing the numbers used in the calculation
    print("Likelihood Ratio Test (Chi-Square Difference Test)")
    print("This test compares the fit of the full model (Model 1) to the constrained model (Model 2).")
    print(f"Chi-Square of Constrained Model (Model 2): {chi_square_model_2:.3f}")
    print(f"Chi-Square of Full Model (Model 1): {chi_square_model_1:.3f}")
    print("---")
    print(f"Final Test Statistic = {chi_square_model_2:.3f} - {chi_square_model_1:.3f} = {test_statistic:.3f}")
    print(f"Degrees of Freedom = {df_diff}")
    print(f"P-value = {p_value}")

solve_task()