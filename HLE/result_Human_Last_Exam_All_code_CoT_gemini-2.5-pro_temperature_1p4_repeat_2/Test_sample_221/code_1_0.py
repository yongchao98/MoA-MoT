import pandas as pd
import numpy as np
from sklearn.preprocessing import scale
from semopy import Model
from scipy.stats import chi2

def solve_sem_anova():
    """
    This function replicates the user's R script to determine the expected
    P-value from an ANOVA comparison of two nested SEM models.
    """
    # Set a seed for reproducibility
    np.random.seed(42)

    # 1. Simulate the dataset as in the R code
    n = 100
    
    # In R, `1:2` supplied to a data.frame of length n > 2 will be recycled.
    # We replicate this behavior with np.tile.
    water_availability = np.tile([1, 2], n // 2)

    data_dict = {
        'water_availability': water_availability,
        'plant_intra_div': np.random.choice([1, 2], size=n, replace=True),
        'plant_inter_div': np.random.choice(np.arange(1, 11), size=n, replace=True)
    }
    
    # The `within` block in R calculates new columns based on existing ones.
    data_dict['insect_richness'] = (data_dict['water_availability'] * 0.01 +
                                    data_dict['plant_intra_div'] * 0.5 +
                                    data_dict['plant_inter_div'] * 1.2 +
                                    np.random.normal(size=n))
    
    # The key step: herbivory is strongly dependent on insect_richness
    data_dict['herbivory'] = (data_dict['insect_richness'] * 3.14 +
                              data_dict['water_availability'] * 0.5 +
                              data_dict['plant_intra_div'] * 0.1 +
                              data_dict['plant_inter_div'] * 0.2 +
                              np.random.normal(size=n))

    example_df = pd.DataFrame(data_dict)

    # Replicate R's scale() function to standardize data (mean=0, sd=1)
    example_scaled = pd.DataFrame(scale(example_df), columns=example_df.columns)

    # 2. Define the two models
    # Model 1: The full model including the path from insect_richness to herbivory
    model_1_spec = """
    herbivory ~ insect_richness + water_availability + plant_intra_div + plant_inter_div
    insect_richness ~ water_availability + plant_intra_div + plant_inter_div
    """

    # Model 2: The nested model, constraining the path to zero
    model_2_spec = """
    herbivory ~ water_availability + plant_intra_div + plant_inter_div
    insect_richness ~ water_availability + plant_intra_div + plant_inter_div
    """

    # 3. Fit the models
    model_1 = Model(model_1_spec)
    fit_1 = model_1.fit(example_scaled)

    model_2 = Model(model_2_spec)
    fit_2 = model_2.fit(example_scaled)
    
    # 4. Perform the Chi-Square Difference Test (what anova() does)
    # The test statistic is the difference in the models' chi-square values
    # In semopy, objective_function_value is proportional to the chi-square value
    lrt_statistic = fit_2.objective_function_value - fit_1.objective_function_value

    # The degrees of freedom for the test is the difference in the models' DoF
    dof_diff = fit_2.dof - fit_1.dof

    # The p-value is calculated from the chi-square distribution
    p_value = chi2.sf(lrt_statistic, dof_diff)

    print("--- Chi-Square Difference Test ---")
    print(f"This test compares the full model (Model 1) with the restricted model (Model 2).")
    print(f"The null hypothesis is that the restriction (removing the path from insect richness to herbivory) does not worsen the model fit.")
    print("\nTest Equation: p_value = P(Chi-Square(df) > statistic)")
    print("\nCalculated values:")
    print(f"Chi-Square Statistic = {lrt_statistic:.4f}")
    print(f"Degrees of Freedom (df) = {dof_diff}")
    print(f"Resulting P-value = {p_value}")

    print("\nConclusion:")
    print("Because the data was generated with a very strong link between insect richness and herbivory,")
    print("the model that omits this link (Model 2) fits the data very poorly compared to the full model.")
    print("This results in a large Chi-Square statistic and a P-value that is practically zero.")

solve_sem_anova()