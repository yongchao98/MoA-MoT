import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.api import anova_lm

# for reproducibility
np.random.seed(123)
n = 100

# Step 1: Simulate the data as described in the R script.
# np.tile repeats the array [1, 2] fifty times to mimic `1:2` with n=100
water_availability = np.tile([1, 2], int(n / 2))
plant_intra_div = np.random.choice([1, 2], n, replace=True)
plant_inter_div = np.random.choice(range(1, 11), n, replace=True)

df = pd.DataFrame({
    'water_availability': water_availability,
    'plant_intra_div': plant_intra_div,
    'plant_inter_div': plant_inter_div
})

# Generate dependent variables based on the specified relationships
df['insect_richness'] = (df['water_availability'] * 0.01 +
                         df['plant_intra_div'] * 0.5 +
                         df['plant_inter_div'] * 1.2 +
                         np.random.normal(size=n))

df['herbivory'] = (df['insect_richness'] * 3.14 +
                   df['water_availability'] * 0.5 +
                   df['plant_intra_div'] * 0.1 +
                   df['plant_inter_div'] * 0.2 +
                   np.random.normal(size=n))

# Step 2: Scale the data, as done in the R script.
# This standardizes the data to have a mean of 0 and a standard deviation of 1.
df_scaled = pd.DataFrame(data=(df - df.mean()) / df.std(), columns=df.columns)

# Step 3: Define and fit the two nested models for 'herbivory'.
# This is the Python equivalent of the models specified for lavaan.

# model_1: The full model, including the path from insect_richness to herbivory
full_model = ols(
    'herbivory ~ insect_richness + water_availability + plant_intra_div + plant_inter_div',
    data=df_scaled
).fit()

# model_2: The reduced model, omitting the path from insect_richness to herbivory
reduced_model = ols(
    'herbivory ~ water_availability + plant_intra_div + plant_inter_div',
    data=df_scaled
).fit()

# Step 4: Perform an ANOVA to compare the two nested models.
# This test determines if the additional variable in the full model (insect_richness)
# provides a statistically significant improvement in fit.
anova_results = anova_lm(reduced_model, full_model)
p_value = anova_results['Pr(>F)'].iloc[1]

# The p-value from this test will be extremely small because the 'insect_richness'
# variable is a very strong predictor of 'herbivory' by design.
# In statistics, a value this small is often represented as 0.0.
print("The anova() function compares the fit of two nested models.")
print("Model 1 (full model) includes the path: insect_richness -> herbivory")
print("Model 2 (reduced model) omits this path.")
print("The data was simulated with a very strong effect (coefficient = 3.14) for this path.")
print("\nTherefore, removing the path will significantly worsen the model fit.")
print("The ANOVA test's null hypothesis is that the reduced model fits as well as the full model.")
print("We expect to strongly reject this null hypothesis, resulting in a very small p-value.")
print(f"\nThe calculated P value from a Python simulation is: {p_value}")
print("This value is effectively 0.")
