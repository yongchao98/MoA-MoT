import pandas as pd
import numpy as np
# The R code uses the 'lavaan' package. The Python equivalent for Structural Equation Modeling is 'semopy'.
# We need to install it if it's not already present: pip install semopy
try:
    from semopy import Model
    from scipy.stats import chi2
except ImportError:
    print("Please install the required libraries: pip install pandas numpy semopy scipy")
    # Exit if libraries are not found
    exit()

# Set a seed for reproducibility of the random data
np.random.seed(42)

# --- 1. Data Generation ---
# This section replicates the data generation process from the R script.
n = 100
# In R, `water_availability = 1:2` is recycled in the data frame creation.
water_availability = np.tile([1, 2], n // 2)

df_raw = pd.DataFrame({
    'plant_intra_div': np.random.randint(1, 3, n),
    'plant_inter_div': np.random.randint(1, 11, n),
    'water_availability': water_availability
})

df_raw['insect_richness'] = (df_raw['water_availability'] * 0.01 +
                             df_raw['plant_intra_div'] * 0.5 +
                             df_raw['plant_inter_div'] * 1.2 +
                             np.random.randn(n))

# The key relationship: herbivory is strongly driven by insect_richness (coefficient = 3.14)
df_raw['herbivory'] = (df_raw['insect_richness'] * 3.14 +
                       df_raw['water_availability'] * 0.5 +
                       df_raw['plant_intra_div'] * 0.1 +
                       df_raw['plant_inter_div'] * 0.2 +
                       np.random.randn(n))

# The R code scales the entire dataset. We will do the same.
df_scaled = (df_raw - df_raw.mean()) / df_raw.std()

# Python-friendly column names
df_scaled.columns = ['plant_intra_div', 'plant_inter_div', 'water_availability', 'insect_richness', 'herbivory']

# --- 2. Model Definition ---
# Model 1: The full model, includes the path from insect richness to herbivory
model_1_spec = """
herbivory ~ insect_richness + water_availability + plant_intra_div + plant_inter_div
insect_richness ~ water_availability + plant_intra_div + plant_inter_div
"""
# Model 2: The nested/restricted model, omits the path from insect richness to herbivory
model_2_spec = """
herbivory ~ water_availability + plant_intra_div + plant_inter_div
insect_richness ~ water_availability + plant_intra_div + plant_inter_div
"""

# --- 3. Model Fitting ---
m1 = Model(model_1_spec)
fit_1 = m1.fit(df_scaled, obj='MLW') # Using Wishart likelihood, similar to lavaan default
m2 = Model(model_2_spec)
fit_2 = m2.fit(df_scaled, obj='MLW')

# --- 4. Model Comparison (Chi-Square Difference Test) ---
# The anova() function in lavaan performs this test. We calculate it here step-by-step.
chi2_1 = fit_1.results['chi2'][0]
dof_1 = fit_1.results['dof'][0]
chi2_2 = fit_2.results['chi2'][0]
dof_2 = fit_2.results['dof'][0]

# The test statistic is the difference in Chi-Square values
chi2_difference = chi2_2 - chi2_1
# The degrees of freedom for the test is the difference in the models' degrees of freedom
dof_difference = dof_2 - dof_1
# Calculate the p-value from the chi-square distribution
p_value = chi2.sf(chi2_difference, dof_difference)

print("--- Chi-Square Difference Test Components ---")
print(f"Model 1 Chi-Square = {chi2_1:.4f}")
print(f"Model 2 Chi-Square = {chi2_2:.4f}")
print(f"Chi-Square Difference = {chi2_difference:.4f}")
print(f"Degrees of Freedom Difference = {dof_difference}")
print("\n--- Final Answer ---")
print(f"The expected P-value is: {p_value}")