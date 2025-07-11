import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.api import lr_test
from scipy.stats import chi2

def solve_sem_puzzle():
    """
    This function replicates the R script's logic in Python to find the
    expected p-value from the anova comparison.
    """
    # Set a seed for reproducibility of the random data
    np.random.seed(42)

    # 1. Recreate the data generation process from the R script
    n = 100
    # In R, `1:2` is recycled to the length of the dataframe.
    water_availability = np.tile([1, 2], n // 2)
    # R's sample(1:2) is equivalent to numpy's randint(1, 3)
    plant_intra_div = np.random.randint(1, 3, size=n)
    # R's sample(1:10) is equivalent to numpy's randint(1, 11)
    plant_inter_div = np.random.randint(1, 11, size=n)

    example_unscaled = pd.DataFrame({
        'water_availability': water_availability,
        'plant_intra_div': plant_intra_div,
        'plant_inter_div': plant_inter_div
    })

    # Generate the dependent variables based on the specified formulas
    insect_richness = (example_unscaled['water_availability'].values * 0.01 +
                       example_unscaled['plant_intra_div'].values * 0.5 +
                       example_unscaled['plant_inter_div'].values * 1.2 +
                       np.random.randn(n))

    herbivory = (insect_richness * 3.14 +
                 example_unscaled['water_availability'].values * 0.5 +
                 example_unscaled['plant_intra_div'].values * 0.1 +
                 example_unscaled['plant_inter_div'].values * 0.2 +
                 np.random.randn(n))

    example_unscaled['insect_richness'] = insect_richness
    example_unscaled['herbivory'] = herbivory

    # Scale all variables, as done in the R script
    scaler = StandardScaler()
    example_scaled = pd.DataFrame(scaler.fit_transform(example_unscaled), columns=example_unscaled.columns)

    # 2. Fit the two models
    # The lavaan anova test on these models is equivalent to a likelihood
    # ratio test between two OLS models predicting 'herbivory'.

    # Model 1 (Full model)
    model_1_formula = 'herbivory ~ insect_richness + water_availability + plant_intra_div + plant_inter_div'
    fit_1 = smf.ols(model_1_formula, data=example_scaled).fit()

    # Model 2 (Reduced model, where `insect_richness` is omitted)
    model_2_formula = 'herbivory ~ water_availability + plant_intra_div + plant_inter_div'
    fit_2 = smf.ols(model_2_formula, data=example_scaled).fit()

    # 3. Perform the Likelihood Ratio Test to compare the models
    lrt_statistic, p_value, df_diff = lr_test(fit_2, fit_1)

    # 4. Print the results, showing the equation for the LRT statistic
    print("The Likelihood Ratio Test (LRT) compares the two nested models.")
    print("The test statistic is calculated as: LRT = 2 * (LogLikelihood_Full - LogLikelihood_Reduced)")
    
    # We use -fit_2.llf because the formula uses the negative of the reduced model's likelihood
    print(f"The final equation is: {lrt_statistic:.4f} = 2 * ({fit_1.llf:.4f} - ({fit_2.llf:.4f}))")
    print("-" * 30)
    print(f"LRT Statistic: {lrt_statistic:.4f}")
    print(f"Degrees of Freedom: {df_diff}")
    print(f"P-value: {p_value}")
    print("-" * 30)
    print("Conclusion: The p-value is extremely small, effectively zero, indicating that the full model is a significantly better fit.")
    print("This is expected because the data was generated with a strong link between 'insect_richness' and 'herbivory'.")


solve_sem_puzzle()