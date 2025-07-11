import pandas as pd
import numpy as np
import semopy
from sklearn.preprocessing import StandardScaler

# Set a random seed for reproducibility, so the result is consistent.
np.random.seed(123)

# --- Step 1: Replicate the R data generation process ---
n = 100
# In R, `1:2` is recycled to fill the 100 rows.
water_availability_vals = np.tile([1, 2], n // 2)

# In R, `sample(1:2, n, replace=T)` samples 1 and 2.
# In R, `sample(1:10, n, replace=T)` samples integers from 1 to 10.
df = pd.DataFrame({
    'water_availability': water_availability_vals,
    'plant_intra_div': np.random.randint(1, 3, n),
    'plant_inter_div': np.random.randint(1, 11, n)
})

# Generate the dependent variables according to the specified formulas.
# The `np.random.randn(n)` call replicates R's `rnorm(n)`.
df['insect_richness'] = (df['water_availability'] * 0.01 +
                         df['plant_intra_div'] * 0.5 +
                         df['plant_inter_div'] * 1.2 +
                         np.random.randn(n))

# Crucially, 'herbivory' is created with a strong dependency on 'insect_richness' (coefficient = 3.14).
df['herbivory'] = (df['insect_richness'] * 3.14 +
                   df['water_availability'] * 0.5 +
                   df['plant_intra_div'] * 0.1 +
                   df['plant_inter_div'] * 0.2 +
                   np.random.randn(n))

# Replicate the `scale()` function from R.
scaler = StandardScaler()
scaled_data = scaler.fit_transform(df)
scaled_df = pd.DataFrame(scaled_data, columns=df.columns)

# --- Step 2: Define and fit the two nested SEM models ---
# Model 1 includes the path from insect_richness to herbivory.
model_1_desc = """
herbivory ~ insect_richness + water_availability + plant_intra_div + plant_inter_div
insect_richness ~ water_availability + plant_intra_div + plant_inter_div
"""

# Model 2 omits this path, effectively constraining its coefficient to 0.
model_2_desc = """
herbivory ~ water_availability + plant_intra_div + plant_inter_div
insect_richness ~ water_availability + plant_intra_div + plant_inter_div
"""

# Fit both models using the semopy library (a Python equivalent to lavaan).
model_1 = semopy.Model(model_1_desc)
fit_1 = model_1.fit(scaled_df)

model_2 = semopy.Model(model_2_desc)
fit_2 = model_2.fit(scaled_df)

# --- Step 3: Compare the models using a Likelihood Ratio Test ---
# The null hypothesis is that the simpler model (Model 2) fits as well as the more complex model (Model 1).
comparison_results = semopy.compare(fit_1, fit_2)

# Extract the key values from the comparison.
chi_sq_m1 = fit_1.calc_gfi() # In semopy, gfi() on a saturated model returns the Chi-square value.
chi_sq_m2 = fit_2.calc_gfi()
chi_sq_diff = comparison_results['Chisq diff'][1]
df_diff = comparison_results['Df diff'][1]
p_value = comparison_results['p-value'][1]

# --- Step 4: Print the results of the comparison test ---
print("Likelihood Ratio Test (Chi-Square Difference Test) Results:")
print(f"The test compares the fit of the full model (Model 1) against the constrained model (Model 2).")
print("-" * 50)
# Using Unicode for delta and chi
print(f"Equation: Δχ² = χ²(Model 2) - χ²(Model 1)")
print(f"Calculation: Δχ² = {chi_sq_m2:.4f} - {chi_sq_m1:.4f} = {chi_sq_diff:.4f}")
print(f"Degrees of Freedom Difference: {df_diff}")
print(f"Resulting P-value: {p_value}")
