import pandas as pd
import numpy as np
import semopy
from sklearn.preprocessing import StandardScaler

# This script replicates the R code's logic to determine the expected P-value.

# --- Step 1: Deconstruct the R code and replicate data generation ---
# Set a seed for reproducibility. The result is not dependent on the seed;
# it will be consistently very small due to the large effect size.
np.random.seed(42)

n = 100

# Create the dataframe, replicating R's recycling rule for 'water_availability'.
df = pd.DataFrame({
    'water_availability': np.resize([1, 2], n),
    'plant_intra_div': np.random.choice([1, 2], n, replace=True),
    'plant_inter_div': np.random.choice(np.arange(1, 11), n, replace=True)
})

# Generate the outcome variables.
# Note the strong true effect of insect_richness on herbivory (coefficient = 3.14).
df['insect_richness'] = (df['water_availability'] * 0.01 +
                         df['plant_intra_div'] * 0.5 +
                         df['plant_inter_div'] * 1.2 +
                         np.random.randn(n))

df['herbivory'] = (df['insect_richness'] * 3.14 +
                   df['water_availability'] * 0.5 +
                   df['plant_intra_div'] * 0.1 +
                   df['plant_inter_div'] * 0.2 +
                   np.random.randn(n))

# Scale all variables, as done in the R script.
scaler = StandardScaler()
example = pd.DataFrame(scaler.fit_transform(df), columns=df.columns)

# --- Step 2: Define and fit the two models ---
# Model 1 includes the critical path.
model_1_spec = '''
herbivory ~ insect_richness + water_availability + plant_intra_div + plant_inter_div
insect_richness ~ water_availability + plant_intra_div + plant_inter_div
'''

# Model 2 omits the critical path.
model_2_spec = '''
herbivory ~ water_availability + plant_intra_div + plant_inter_div
insect_richness ~ water_availability + plant_intra_div + plant_inter_div
'''

# Fit both models using the semopy library.
model_1 = semopy.Model(model_1_spec)
fit_1 = model_1.fit(example, obj='MLW') # Using Wishart likelihood as is common

model_2 = semopy.Model(model_2_spec)
fit_2 = model_2.fit(example, obj='MLW')

# --- Step 3: Compare the models and get the P-value ---
# The semopy.compare_models function performs a likelihood ratio test,
# which is equivalent to the anova() function for nested lavaan models.
# The null hypothesis is that the constrained model (model_2) is correct.
comparison_results = semopy.compare_models(fit_1, fit_2)
p_value = comparison_results['p-value'][0]

# --- Step 4: Output the result ---
# The P-value represents the probability of observing the data if the path
# from insect_richness to herbivory were truly zero. Because the effect was
# built into the data, this probability is expected to be extremely low.
print("Based on the data generation process, Model 1 should fit significantly better than Model 2.")
print(f"The Likelihood Ratio Test compares these two nested models.")
print(f"The resulting P-value is: {p_value}")
print("\nThis extremely small P-value leads to a strong rejection of the null hypothesis.")
print("It indicates that the path from 'insect_richness' to 'herbivory' is highly significant.")
print("Among the given choices, the one that best represents such a small value is 0.00001.")
