import pandas as pd
import numpy as np
import semopy
from sklearn.preprocessing import StandardScaler
import warnings

# Suppress convergence warnings from semopy for this example
warnings.filterwarnings("ignore", category=semopy.estimator.ConvergenceWarning)

# 1. Set a seed for reproducibility
np.random.seed(42)

# 2. Generate Data (replicating the R script's logic)
n = 100
# R's recycling of `1:2` results in alternating 1s and 2s for water_availability
df = pd.DataFrame({
    'water_availability': np.tile([1, 2], n // 2),
    # R's sample(1:2) is inclusive, so randint(1,3) is the Python equivalent
    'plant_intra_div': np.random.randint(1, 3, n),
    # R's sample(1:10) is inclusive, so randint(1,11) is the Python equivalent
    'plant_inter_div': np.random.randint(1, 11, n)
})

# Calculate dependent variables based on the formulas in the problem
df['insect_richness'] = (df['water_availability'] * 0.01 +
                         df['plant_intra_div'] * 0.5 +
                         df['plant_inter_div'] * 1.2 +
                         np.random.randn(n))

df['herbivory'] = (df['insect_richness'] * 3.14 +
                   df['water_availability'] * 0.5 +
                   df['plant_intra_div'] * 0.1 +
                   df['plant_inter_div'] * 0.2 +
                   np.random.randn(n))

# 3. Scale the final DataFrame, as in the R script
scaler = StandardScaler()
example_df = pd.DataFrame(scaler.fit_transform(df), columns=df.columns)

# 4. Define the two models using semopy syntax
# Model 1 includes the path from insect_richness to herbivory
model_1_desc = '''
herbivory ~ insect_richness + water_availability + plant_intra_div + plant_inter_div
insect_richness ~ water_availability + plant_intra_div + plant_inter_div
'''

# Model 2 omits the path from insect_richness to herbivory
model_2_desc = '''
herbivory ~ water_availability + plant_intra_div + plant_inter_div
insect_richness ~ water_availability + plant_intra_div + plant_inter_div
'''

# 5. Fit both models to the data
fit_1 = semopy.Model(model_1_desc)
res_1 = fit_1.fit(example_df)

fit_2 = semopy.Model(model_2_desc)
res_2 = fit_2.fit(example_df)

# 6. Compare the models using a Chi-Squared Difference Test
# The result is stored in a pandas DataFrame
comparison_results = semopy.compare(fit_1, fit_2)
result_row = comparison_results.loc['anova_1']

# 7. Print the results of the comparison
print("Model comparison results (Chi-Squared Difference Test):")
# The "final equation" in this context is the test statistic and its p-value.
# Here are the numbers that make up that result.
print(f"Chi-Squared Difference: {result_row['Chisq-diff']:.4f}")
print(f"Degrees of Freedom Difference: {int(result_row['df-diff'])}")
print(f"P-value: {result_row['p-value']:.10f}")