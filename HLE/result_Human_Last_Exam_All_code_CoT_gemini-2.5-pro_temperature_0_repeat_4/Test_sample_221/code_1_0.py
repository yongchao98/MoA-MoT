import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf

# Set a seed for reproducibility
np.random.seed(42)

# 1. Simulate data similar to the R script
n = 100
# In the R code, water_availability is 1:2, but this will be perfectly collinear with the intercept
# after scaling if not repeated. We'll simulate it as a random binary variable.
example_df = pd.DataFrame({
    'water_availability': np.random.randint(1, 3, n),
    'plant_intra_div': np.random.randint(1, 3, n),
    'plant_inter_div': np.random.randint(1, 11, n)
})

# Create the dependent variables based on the formulas
example_df['insect_richness'] = (example_df['water_availability'] * 0.01 +
                                 example_df['plant_intra_div'] * 0.5 +
                                 example_df['plant_inter_div'] * 1.2 +
                                 np.random.randn(n))

example_df['herbivory'] = (example_df['insect_richness'] * 3.14 +
                           example_df['water_availability'] * 0.5 +
                           example_df['plant_intra_div'] * 0.1 +
                           example_df['plant_inter_div'] * 0.2 +
                           np.random.randn(n))

# Scale the data (mean=0, std=1) as in the R script
example_df_scaled = pd.DataFrame(
    data= (example_df - example_df.mean()) / example_df.std(),
    columns=example_df.columns
)


# 2. Define and fit two models corresponding to model_1 and model_2
# Model 1: Full model including the insect_richness -> herbivory path
model_1_formula = 'herbivory ~ insect_richness + water_availability + plant_intra_div + plant_inter_div'
fit_1 = smf.ols(model_1_formula, data=example_df_scaled).fit()

# Model 2: Reduced model, omitting the insect_richness -> herbivory path
model_2_formula = 'herbivory ~ water_availability + plant_intra_div + plant_inter_div'
fit_2 = smf.ols(model_2_formula, data=example_df_scaled).fit()

# 3. Compare the two nested models using ANOVA
# This tests the null hypothesis that the simpler model (fit_2) is sufficient.
anova_results = sm.stats.anova_lm(fit_2, fit_1)

# The P-value is in the 'Pr(>F)' column of the second row
p_value = anova_results.iloc[1]['Pr(>F)']

# 4. Print the result
print("The analysis compares two models:")
print("Model 1: herbivory ~ insect_richness + water_availability + plant_intra_div + plant_inter_div")
print("Model 2: herbivory ~ water_availability + plant_intra_div + plant_inter_div")
print("\nThe data was generated with a strong effect of 'insect_richness' on 'herbivory' (coefficient = 3.14).")
print("The ANOVA test checks if removing this effect (as in Model 2) significantly worsens the model fit.")
print(f"\nBecause the effect is very strong, we expect a highly significant result (a very small P-value).")
print("\n--- ANOVA Results ---")
print(anova_results)
print("\n--- Conclusion ---")
print(f"The P-value from the comparison is: {p_value}")
print("This value is extremely small, confirming that Model 1 is a significantly better fit.")
print("Among the given choices, 0.00001 is the most plausible answer.")
