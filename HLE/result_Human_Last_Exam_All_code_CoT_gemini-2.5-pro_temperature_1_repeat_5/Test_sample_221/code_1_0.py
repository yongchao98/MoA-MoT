import numpy as np
import pandas as pd
from sklearn.preprocessing import scale
import statsmodels.api as sm
from statsmodels.stats.api import anova_lm

# Set a seed for reproducibility
np.random.seed(42)

# 1. Simulate the dataset as described in the R code
n = 100
# In R, 1:2 is recycled to fill 100 spots: 1, 2, 1, 2, ...
water_availability = np.tile([1, 2], n // 2)

df = pd.DataFrame({
    'water_availability': water_availability,
    'plant_intra_div': np.random.randint(1, 3, n),  # R's sample(1:2, n, replace=T)
    'plant_inter_div': np.random.randint(1, 11, n) # R's sample(1:10, n, replace=T)
})

# Create the dependent variables based on the specified relationships
df['insect_richness'] = (df['water_availability'] * 0.01 +
                         df['plant_intra_div'] * 0.5 +
                         df['plant_inter_div'] * 1.2 +
                         np.random.randn(n)) # rnorm(n)

df['herbivory'] = (df['insect_richness'] * 3.14 + # The key strong relationship
                   df['water_availability'] * 0.5 +
                   df['plant_intra_div'] * 0.1 +
                   df['plant_inter_div'] * 0.2 +
                   np.random.randn(n))

# Scale the data, as in the R code
scaled_df = pd.DataFrame(scale(df), columns=df.columns)

# 2. Define and fit models analogous to the SEM models
# The core question is whether 'insect_richness' significantly predicts 'herbivory'

# Model 1 (Full model)
Y = scaled_df['herbivory']
X_full = scaled_df[['insect_richness', 'water_availability', 'plant_intra_div', 'plant_inter_div']]
X_full = sm.add_constant(X_full) # Add intercept
model_1_fit = sm.OLS(Y, X_full).fit()

# Model 2 (Reduced model, without insect_richness)
X_reduced = scaled_df[['water_availability', 'plant_intra_div', 'plant_inter_div']]
X_reduced = sm.add_constant(X_reduced) # Add intercept
model_2_fit = sm.OLS(Y, X_reduced).fit()

# 3. Compare the nested models to get a P-value
# anova_lm performs an F-test, which is the OLS equivalent of the chi-square difference test
anova_results = anova_lm(model_2_fit, model_1_fit)
p_value = anova_results['Pr(>F)'][1]

print("The P-value from comparing the models is:")
print(p_value)