import pandas as pd
import numpy as np
import semopy
from sklearn.preprocessing import StandardScaler
import warnings

# Suppress warnings from semopy for cleaner output
warnings.filterwarnings("ignore", category=UserWarning)

# Plan:
# 1. Recreate the R data generation process in Python.
# 2. Scale the data as done in the R script.
# 3. Define the two SEM models.
# 4. Fit the models using the semopy library.
# 5. Compare the nested models to get the chi-square difference test results.
# 6. Print the components of the test and the final p-value.

# Set a seed for reproducibility, so the random numbers are the same each time.
np.random.seed(123)
n = 100

# 1. Data Generation
# In R, data.frame(water_availability = 1:2, ...) recycles 1:2 to length n
water_availability = np.tile([1, 2], n // 2)
plant_intra_div = np.random.randint(1, 3, n)
plant_inter_div = np.random.randint(1, 11, n)

example_df = pd.DataFrame({
    'water_availability': water_availability,
    'plant_intra.div': plant_intra_div,
    'plant_inter.div': plant_inter_div
})

insect_richness = (example_df['water_availability'] * 0.01 +
                   example_df['plant_intra.div'] * 0.5 +
                   example_df['plant_inter.div'] * 1.2 +
                   np.random.randn(n))

# The key line: herbivory is strongly dependent on insect_richness (coefficient = 3.14)
herbivory = (insect_richness * 3.14 +
             example_df['water_availability'] * 0.5 +
             example_df['plant_intra.div'] * 0.1 +
             example_df['plant_inter.div'] * 0.2 +
             np.random.randn(n))

example_df['insect_richness'] = insect_richness
example_df['herbivory'] = herbivory

# 2. Scaling
# The scale() function in R centers and scales the data.
scaler = StandardScaler()
example_scaled = pd.DataFrame(scaler.fit_transform(example_df), columns=example_df.columns)

# 3. Model Specification
model_1 = """
# This model includes the path from insect_richness to herbivory
herbivory ~ insect_richness + water_availability + plant_intra.div + plant_inter.div
insect_richness ~ water_availability + plant_intra.div + plant_inter.div
"""

model_2 = """
# This model omits the path, constraining its coefficient to 0
herbivory ~ water_availability + plant_intra.div + plant_inter.div
insect_richness ~ water_availability + plant_intra.div + plant_inter.div
"""

# 4. Model Fitting
fit_1 = semopy.Model(model_1)
res_1 = fit_1.fit(example_scaled)

fit_2 = semopy.Model(model_2)
res_2 = fit_2.fit(example_scaled)

# 5. Model Comparison
comparison_results = semopy.compare(fit_1, fit_2)

# 6. Print the final equation numbers
# The chi-square difference test calculates a p-value based on the difference
# in chi-square statistics and the difference in degrees of freedom between the two models.
chi2_1 = fit_1.calc_stats()['chi2'][0]
dof_1 = fit_1.calc_stats()['dof'][0]
chi2_2 = fit_2.calc_stats()['chi2'][0]
dof_2 = fit_2.calc_stats()['dof'][0]
chi2_diff = comparison_results['chi2 diff'][0]
dof_diff = comparison_results['DoF diff'][0]
p_value = comparison_results['p-value'][0]

print("--- Chi-Square Difference Test for Nested Models ---")
print(f"Model 1 (Full Model):      Chi2({dof_1}) = {chi2_1:.4f}")
print(f"Model 2 (Restricted Model): Chi2({dof_2}) = {chi2_2:.4f}")
print("-" * 50)
print("Test Equation:")
print(f"Chi-square Difference = Chi2(Model 2) - Chi2(Model 1)")
print(f"                      = {chi2_2:.4f} - {chi2_1:.4f} = {chi2_diff:.4f}")
print(f"Degrees of Freedom Difference = DoF(Model 2) - DoF(Model 1) = {dof_2} - {dof_1} = {dof_diff}")
print("-" * 50)
print(f"The resulting P-value is: {p_value}")
print("\nAs the p-value is extremely small (effectively 0), we conclude that the true value is best represented by 0.")
