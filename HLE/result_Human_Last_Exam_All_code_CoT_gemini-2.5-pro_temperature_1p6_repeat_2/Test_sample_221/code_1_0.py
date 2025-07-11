import numpy as np
import pandas as pd
from sklearn.preprocessing import scale
import statsmodels.api as sm
from scipy.stats import chi2

# Set a seed for reproducibility of random numbers
np.random.seed(42)

# 1. Simulate the dataset as described in the R code
n = 100
df = pd.DataFrame({
    # Use np.tile to replicate the 1:2 pattern for n=100
    'water_availability': np.tile([1, 2], n // 2),
    'plant_intra_div': np.random.randint(1, 3, n),
    'plant_inter_div': np.random.randint(1, 11, n)
})

# Generate the dependent variables based on the formulas
df['insect_richness'] = (df['water_availability'] * 0.01 +
                         df['plant_intra_div'] * 0.5 +
                         df['plant_inter_div'] * 1.2 +
                         np.random.normal(size=n))

true_coefficient = 3.14
df['herbivory'] = (df['insect_richness'] * true_coefficient +
                   df['water_availability'] * 0.5 +
                   df['plant_intra_div'] * 0.1 +
                   df['plant_inter_div'] * 0.2 +
                   np.random.normal(size=n))

# Scale all variables as done in the R code
df_scaled = pd.DataFrame(scale(df), columns=df.columns)

# Define the models for the likelihood ratio test
Y = df_scaled['herbivory']

# 2. Fit Model 1 (Full Model), which includes the disputed path
X1 = df_scaled[['insect_richness', 'water_availability', 'plant_intra_div', 'plant_inter_div']]
X1 = sm.add_constant(X1)
model1 = sm.OLS(Y, X1).fit()
loglik1 = model1.llf

# 3. Fit Model 2 (Reduced Model), which excludes the disputed path
X2 = df_scaled[['water_availability', 'plant_intra_div', 'plant_inter_div']]
X2 = sm.add_constant(X2)
model2 = sm.OLS(Y, X2).fit()
loglik2 = model2.llf

# 4. Perform the Likelihood Ratio Test
# The number of parameters differs by 1 (the 'insect_richness' coefficient)
df_diff = model1.df_model - model2.df_model
# The test statistic is 2 * (log-likelihood of full model - log-likelihood of reduced model)
lr_statistic = 2 * (loglik1 - loglik2)
# The p-value is calculated from the chi-square distribution
p_value = chi2.sf(lr_statistic, df_diff)

print(f"The 'true' coefficient for the 'insect_richness -> herbivory' path is: {true_coefficient}")
print(f"The 'anova' test checks if this coefficient is significantly different from 0.")
print(f"Likelihood Ratio Test Statistic (Chi-Square): {lr_statistic:.4f}")
print(f"Degrees of Freedom: {df_diff:.0f}")
print(f"Resulting P-value: {p_value}")

# Since the p-value is extremely small, it is effectively 0.
print("\nThe expected P-value is effectively 0.")