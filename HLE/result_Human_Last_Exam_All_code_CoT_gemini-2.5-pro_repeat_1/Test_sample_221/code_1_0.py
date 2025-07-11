import pandas as pd
import numpy as np
from statsmodels.formula.api import ols
from scipy.stats import chi2

# Set a seed for reproducibility to get a consistent result.
np.random.seed(42)

# Step 1: Replicate the data generation process from the R script.
n = 100
# In R, `1:2` gets recycled to match the length of n.
water_availability = np.tile(np.array([1, 2]), n // 2)
df = pd.DataFrame({
    'water_availability': water_availability,
    'plant_intra_div': np.random.randint(1, 3, n),
    'plant_inter_div': np.random.randint(1, 11, n)
})

# Generate the dependent variables according to the specified formulas.
# Note the strong, direct effect of insect_richness on herbivory.
df['insect_richness'] = (df['water_availability'] * 0.01 +
                         df['plant_intra_div'] * 0.5 +
                         df['plant_inter_div'] * 1.2 +
                         np.random.randn(n))

df['herbivory'] = (df['insect_richness'] * 3.14 +  # This is the crucial path
                   df['water_availability'] * 0.5 +
                   df['plant_intra_div'] * 0.1 +
                   df['plant_inter_div'] * 0.2 +
                   np.random.randn(n))

# The R code scales the data, which standardizes variables to mean=0, sd=1.
# This does not change the underlying relationships, only their scale.
df_scaled = (df - df.mean()) / df.std()

# Step 2: Fit two linear models analogous to the two SEMs being compared.
# The anova in R for these SEMs tests the significance of the omitted path,
# which is equivalent to an ANOVA comparing a full vs. reduced linear model.

# Model 1 (full model), analogous to R `model_1`
model_1_formula = 'herbivory ~ insect_richness + water_availability + plant_intra_div + plant_inter_div'
fit_1 = ols(model_1_formula, data=df_scaled).fit()

# Model 2 (reduced model), analogous to R `model_2`
model_2_formula = 'herbivory ~ water_availability + plant_intra_div + plant_inter_div'
fit_2 = ols(model_2_formula, data=df_scaled).fit()

# Step 3: Perform a Likelihood Ratio (LR) test to compare the models.
# This is what `lavaan::anova` does. The test statistic is distributed as Chi-squared.
loglik_full = fit_1.llf
loglik_reduced = fit_2.llf
# Degrees of freedom for the test is the difference in the number of estimated parameters.
df_difference = fit_1.df_model - fit_2.df_model

lr_statistic = 2 * (loglik_full - loglik_reduced)
p_value = chi2.sf(lr_statistic, df_difference)

# Step 4: Print the results of the comparison.
print("--- Model Comparison using Likelihood Ratio Test ---")
print("This test is analogous to the `anova(fit_1, fit_2)` in the R script.")
print(f"H0: The model without the 'insect_richness -> herbivory' path fits as well as the full model.")
print(f"H1: The full model fits significantly better.")
print("\nTest Results:")
print(f"Chi-squared statistic = {lr_statistic:.4f}")
print(f"Degrees of freedom = {df_difference:.0f}")
print(f"P-value = {p_value}")
print("\nConclusion: Since the P-value is extremely small (effectively 0), we reject the null hypothesis.")
print("The expected P-value from the R function is 0.")
