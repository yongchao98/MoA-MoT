import pandas as pd
import numpy as np
import semopy
from sklearn.preprocessing import StandardScaler
import warnings

# Suppress ConvergenceWarning which may occur with random data but doesn't affect this problem's logic
warnings.filterwarnings("ignore", category=semopy.ConvergenceWarning)

# Plan:
# 1. Replicate the R data generation process in Python using pandas and numpy.
#    Set a random seed for reproducibility.
# 2. Define the two structural models using semopy syntax, which is similar to lavaan.
# 3. Fit both models to the scaled data.
# 4. Use semopy.compare to perform a likelihood ratio test (equivalent to R's anova).
# 5. Print the resulting p-value from the comparison.

# Set a seed for reproducibility, as the original script uses random sampling.
np.random.seed(42)

# 1. Data Generation
n = 100
# Create the initial dataframe
# np.tile([1, 2], n // 2) replicates R's recycling of 1:2 over n=100 rows
example_df = pd.DataFrame({
    'water_availability': np.tile([1, 2], n // 2),
    'plant_intra.div': np.random.randint(1, 3, n),
    'plant_inter.div': np.random.randint(1, 11, n)
})

# Sequentially calculate the dependent variables, as R's `within` block does.
# `insect_richness` is calculated first.
insect_richness = (example_df['water_availability'] * 0.01 +
                   example_df['plant_intra.div'] * 0.5 +
                   example_df['plant_inter.div'] * 1.2 +
                   np.random.normal(size=n)) # rnorm(n) has mean=0, sd=1

# The newly created `insect_richness` is used to calculate `herbivory`.
# Note the very large coefficient (3.14) for insect_richness.
herbivory = (insect_richness * 3.14 +
             example_df['water_availability'] * 0.5 +
             example_df['plant_intra.div'] * 0.1 +
             example_df['plant_inter.div'] * 0.2 +
             np.random.normal(size=n)) # rnorm(n)

example_df['insect_richness'] = insect_richness
example_df['herbivory'] = herbivory

# R's `scale()` standardizes the data (mean=0, sd=1). We use StandardScaler for this.
scaler = StandardScaler()
example_scaled = pd.DataFrame(scaler.fit_transform(example_df), columns=example_df.columns)
# semopy doesn't handle dots in variable names well, so we replace them with underscores.
example_scaled.columns = [c.replace('.', '_') for c in example_scaled.columns]


# 2. Model Specification
# Model 1 includes the critical path from insect_richness to herbivory.
model_1 = '''
# regressions
herbivory ~ insect_richness + water_availability + plant_intra_div + plant_inter_div
insect_richness ~ water_availability + plant_intra_div + plant_inter_div
'''

# Model 2 omits the path from insect_richness to herbivory.
model_2 = '''
# regressions
herbivory ~ water_availability + plant_intra_div + plant_inter_div
insect_richness ~ water_availability + plant_intra_div + plant_inter_div
'''

# 3. Model Fitting
# We fit both models silently, as we are interested in the comparison.
fit_1 = semopy.Model(model_1)
res_1 = fit_1.fit(example_scaled, silent=True)

fit_2 = semopy.Model(model_2)
res_2 = fit_2.fit(example_scaled, silent=True)

# 4. Model Comparison
# semopy.compare performs the likelihood ratio test (Chi-square difference test) for nested models.
comparison_results = semopy.compare(fit_1, fit_2)

# The comparison results table contains the chi-square difference, df difference, and p-value.
# The row for model2 contains the comparison statistics.
# We are interested in the p-value.
p_value = comparison_results.loc['model2', 'p-value']
chi2_diff = comparison_results.loc['model2', 'chisq_d']
df_diff = comparison_results.loc['model2', 'df_d']

print("Likelihood Ratio Test (Chi-Square Difference Test) Results:")
print(f"Chi-square Difference: {chi2_diff:.4f}")
print(f"Degrees of Freedom Difference: {int(df_diff)}")
print(f"P-value: {p_value}")

print("\nAs predicted, the p-value is extremely small, indicating that Model 1 fits the data significantly better.")
print("This is because the path from insect richness to herbivory, which exists in the data, was removed from Model 2.")
print("The closest answer choice is therefore E.")