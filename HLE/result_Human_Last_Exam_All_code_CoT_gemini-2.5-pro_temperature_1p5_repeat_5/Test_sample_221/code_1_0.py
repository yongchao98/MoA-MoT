import pandas as pd
import numpy as np
import semopy

# Set a seed for reproducibility
np.random.seed(42)

# 1. Simulate the data as described in the R code
n = 100
# In R, 1:2 is recycled. The equivalent in numpy is to tile the array.
water_availability_vals = np.tile([1, 2], int(np.ceil(n / 2)))[:n]

df = pd.DataFrame({
    'water_availability': water_availability_vals,
    'plant_intra_div': np.random.choice([1, 2], n, replace=True),
    'plant_inter_div': np.random.randint(1, 11, size=n)
})

# Generate dependent variables based on the formulas
df['insect_richness'] = (df['water_availability'] * 0.01 + 
                         df['plant_intra_div'] * 0.5 + 
                         df['plant_inter_div'] * 1.2 + 
                         np.random.normal(size=n))

df['herbivory'] = (df['insect_richness'] * 3.14 + 
                   df['water_availability'] * 0.5 + 
                   df['plant_intra_div'] * 0.1 + 
                   df['plant_inter_div'] * 0.2 + 
                   np.random.normal(size=n))

# Scale the data
example_scaled = (df - df.mean()) / df.std()

# 2. Define the two models in semopy syntax
model_1 = """
herbivory ~ insect_richness + water_availability + plant_intra_div + plant_inter_div
insect_richness ~ water_availability + plant_intra_div + plant_inter_div
"""

model_2 = """
herbivory ~ water_availability + plant_intra_div + plant_inter_div
insect_richness ~ water_availability + plant_intra_div + plant_inter_div
"""

# 3. Fit the models
fit_1 = semopy.Model(model_1)
res_1 = fit_1.fit(example_scaled)

fit_2 = semopy.Model(model_2)
res_2 = fit_2.fit(example_scaled)

# 4. Compare the models using a chi-square difference test
comparison = semopy.compare(fit_1, fit_2)

# Extract values for the "equation"
chi2_m1 = comparison.loc['model1', 'chi2']
chi2_m2 = comparison.loc['model2', 'chi2']
chi2_diff = chi2_m2 - chi2_m1

dof_m1 = comparison.loc['model1', 'dof']
dof_m2 = comparison.loc['model2', 'dof']
dof_diff = dof_m2 - dof_m1

pval = comparison.loc['model2', 'p-value']

print("Model Comparison (Chi-Square Difference Test)")
print("-" * 45)
print(f"Chi-Square (Full Model 1): {chi2_m1:.4f}")
print(f"Chi-Square (Restricted Model 2): {chi2_m2:.4f}")
print(f"Degrees of Freedom (Full Model 1): {int(dof_m1)}")
print(f"Degrees of Freedom (Restricted Model 2): {int(dof_m2)}")
print("\n--- Test Equation ---")
print(f"Chi-Square Difference = Chi2_Model2 - Chi2_Model1")
print(f"                      = {chi2_m2:.4f} - {chi2_m1:.4f} = {chi2_diff:.4f}")
print(f"DF Difference = {int(dof_diff)}")
print("---------------------\n")
print(f"Resulting P-value: {pval}")
