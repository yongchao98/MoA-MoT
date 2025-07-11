import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.stats.anova import anova_lm

# Set a random seed for reproducible results
np.random.seed(42)

# Step 1: Simulate the data following the logic of the R script
n = 100

# R's '1:2' gets recycled, so we get 1, 2, 1, 2, ...
water_availability = np.tile([1, 2], n // 2)

# sample(1:2, n, replace=T) -> random integers 1 or 2
plant_intra_div = np.random.randint(1, 3, size=n)
# sample(1:10, n, replace=T) -> random integers from 1 to 10
plant_inter_div = np.random.randint(1, 11, size=n)

# Create an initial dataframe
example = pd.DataFrame({
    'water_availability': water_availability,
    'plant_intra_div': plant_intra_div,
    'plant_inter_div': plant_inter_div
})

# Step 2: Calculate the dependent variables based on the formulas
# Generate random noise equivalent to R's rnorm(n)
rnorm_1 = np.random.randn(n)
example['insect_richness'] = (example['water_availability'] * 0.01 +
                            example['plant_intra_div'] * 0.5 +
                            example['plant_inter_div'] * 1.2 + rnorm_1)

rnorm_2 = np.random.randn(n)
example['herbivory'] = (example['insect_richness'] * 3.14 +
                      example['water_availability'] * 0.5 +
                      example['plant_intra_div'] * 0.1 +
                      example['plant_inter_div'] * 0.2 + rnorm_2)

# Step 3: Scale the data, just as R's scale() function does
# It standardizes to mean=0, std=1, using N-1 in the denominator for std
example_scaled = (example - example.mean()) / example.std(ddof=1)

# Step 4: Define and fit the two models to be compared
# Model 2 (Reduced): herbivory ~ water_availability + plant_intra_div + plant_inter_div
X_reduced = example_scaled[['water_availability', 'plant_intra_div', 'plant_inter_div']]
X_reduced = sm.add_constant(X_reduced) # Add intercept
y = example_scaled['herbivory']
fit_2 = sm.OLS(y, X_reduced).fit()

# Model 1 (Full): herbivory ~ insect_richness + water_availability + plant_intra_div + plant_inter_div
X_full = example_scaled[['insect_richness', 'water_availability', 'plant_intra_div', 'plant_inter_div']]
X_full = sm.add_constant(X_full) # Add intercept
fit_1 = sm.OLS(y, X_full).fit()

# Step 5: Perform ANOVA to compare the models
# This tests if adding 'insect_richness' (going from fit_2 to fit_1) significantly improves the model
anova_results = anova_lm(fit_2, fit_1)
p_value = anova_results['Pr(>F)'][1]

# The 'final equation' in this context is the result of the statistical test
f_statistic = anova_results['F'][1]
df_diff = anova_results['df_diff'][1]
ss_diff = anova_results['ss_diff'][1]

print("ANOVA Model Comparison:")
print("The test compares the restricted model (Model 2) against the full model (Model 1).")
print(f"The addition of 'insect_richness' to the model reduces the residual sum of squares by {ss_diff:.4f}.")
print(f"This comparison is based on an F-statistic with {df_diff:.0f} and {fit_1.df_resid:.0f} degrees of freedom.")
print(f"F-Statistic: {f_statistic:.4f}")
print(f"P-value: {p_value}")
