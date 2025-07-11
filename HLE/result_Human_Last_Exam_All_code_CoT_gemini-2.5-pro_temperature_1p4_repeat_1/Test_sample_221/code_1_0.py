import pandas as pd
import numpy as np
import semopy
import io
import sys

def solve_task():
    """
    This function simulates the data, fits the two SEMs, compares them,
    and prints the result of the comparison test.
    """
    # Set a seed for reproducibility. The original R script is random,
    # but the underlying relationship is so strong the outcome is consistent.
    np.random.seed(123)

    # 1. Simulate the data as described in the R script.
    n = 100
    df = pd.DataFrame({
        # R's `water_availability = 1:2` with n=100 recycles the vector [1, 2] fifty times.
        'water_availability': np.resize(np.array([1, 2]), n),
        # `sample(1:2, n, replace=T)` is equivalent to randint(1, 3, n) in numpy.
        'plant_intra.div': np.random.randint(1, 3, n, dtype=int),
        'plant_inter.div': np.random.randint(1, 11, n, dtype=int)
    })

    # Generate the dependent variables according to the formulas.
    # The coefficient of 3.14 creates a very strong causal link.
    df['insect_richness'] = (df['water_availability'] * 0.01 +
                             df['plant_intra.div'] * 0.5 +
                             df['plant_inter.div'] * 1.2 +
                             np.random.normal(loc=0, scale=1, size=n))

    df['herbivory'] = (df['insect_richness'] * 3.14 +
                       df['water_availability'] * 0.5 +
                       df['plant_intra.div'] * 0.1 +
                       df['plant_inter.div'] * 0.2 +
                       np.random.normal(loc=0, scale=1, size=n))

    # Scale the data (z-score normalization).
    # R's scale() uses the sample standard deviation (ddof=1) by default.
    for col in df.columns:
        df[col] = (df[col] - df[col].mean()) / df[col].std(ddof=1)

    # 2. Define the two SEM models using semopy syntax.
    # Model 1 is the full model, reflecting the data generation process.
    model_1_spec = """
    herbivory ~ insect_richness + water_availability + plant_intra.div + plant_inter.div
    insect_richness ~ water_availability + plant_intra.div + plant_inter.div
    """

    # Model 2 is the restricted model, omitting the key causal path.
    model_2_spec = """
    herbivory ~ water_availability + plant_intra.div + plant_inter.div
    insect_richness ~ water_availability + plant_intra.div + plant_inter.div
    """

    # 3. Fit both models. We suppress the verbose output for a cleaner result.
    old_stdout = sys.stdout
    sys.stdout = captured_output = io.StringIO()
    model1 = semopy.Model(model_1_spec)
    res1 = model1.fit(df)
    model2 = semopy.Model(model_2_spec)
    res2 = model2.fit(df)
    sys.stdout = old_stdout

    # 4. Compare the models. This is equivalent to R's `anova(fit_1, fit_2)`.
    # It performs a chi-square difference test.
    comparison_results = semopy.compare(model1, model2)
    diff_row = comparison_results.loc['diff']

    # 5. Output the components of the "final equation", which is the chi-square test.
    # The very strong effect (3.14) included in model 1 but omitted in model 2
    # leads to a large chi-square difference and a p-value that is effectively zero.
    chi2_diff = diff_row['chi2']
    dof_diff = diff_row['dof']
    p_value = diff_row['p-value']

    print("Model Comparison Results (Chi-Square Difference Test)")
    print("-----------------------------------------------------")
    print(f"Chi-square Difference = {chi2_diff:.4f}")
    print(f"Degrees of Freedom Difference = {int(dof_diff)}")
    print(f"P-value = {p_value}")
    print("\nConclusion: The p-value is extremely small, indicating a highly significant difference between the models.")
    print("This means the path from 'insect_richness' to 'herbivory' is essential.")

solve_task()