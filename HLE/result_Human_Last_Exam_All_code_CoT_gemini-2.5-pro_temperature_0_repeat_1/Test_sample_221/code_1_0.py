import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler
from scipy.stats import chi2

# for reproducibility
np.random.seed(42)

n = 100

# 1. Simulate data with a known strong relationship
df = pd.DataFrame({
    'water_availability': np.tile([1, 2], n // 2),
    'plant_intra_div': np.random.randint(1, 3, n),
    'plant_inter_div': np.random.randint(1, 11, n)
})

df['insect_richness'] = (df['water_availability'] * 0.01 +
                         df['plant_intra_div'] * 0.5 +
                         df['plant_inter_div'] * 1.2 +
                         np.random.randn(n))

# The key relationship: herbivory is strongly predicted by insect_richness
df['herbivory'] = (df['insect_richness'] * 3.14 +
                   df['water_availability'] * 0.5 +
                   df['plant_intra_div'] * 0.1 +
                   df['plant_inter_div'] * 0.2 +
                   np.random.randn(n))

# Scale data as in the R script
scaler = StandardScaler()
scaled_data = scaler.fit_transform(df)
scaled_df = pd.DataFrame(scaled_data, columns=df.columns)

# 2. Define and fit the two models
# Model 1 (Full model)
Y = scaled_df['herbivory']
X1 = scaled_df[['insect_richness', 'water_availability', 'plant_intra_div', 'plant_inter_div']]
model_1 = sm.OLS(Y, sm.add_constant(X1)).fit()

# Model 2 (Nested/Reduced model, omits 'insect_richness')
X2 = scaled_df[['water_availability', 'plant_intra_div', 'plant_inter_div']]
model_2 = sm.OLS(Y, sm.add_constant(X2)).fit()

# 3. Perform the Likelihood Ratio Test
# The test statistic is 2 * (log-likelihood of full model - log-likelihood of reduced model)
# The degrees of freedom is the difference in the number of parameters (1 in this case)
log_likelihood_1 = model_1.llf
log_likelihood_2 = model_2.llf
lr_statistic = 2 * (log_likelihood_1 - log_likelihood_2)
df_diff = model_1.df_model - model_2.df_model # Difference in degrees of freedom
p_value = chi2.sf(lr_statistic, df_diff) # sf is the survival function (1 - cdf)

# 4. Print the result
# The anova function in R would output a table, but the key value is the p-value.
# We will print the components of the test and the final p-value.
print(f"Log-Likelihood of Model 1 (full): {log_likelihood_1:.4f}")
print(f"Log-Likelihood of Model 2 (nested): {log_likelihood_2:.4f}")
print(f"Chi-Square Difference Statistic: {lr_statistic:.4f}")
print(f"Difference in Degrees of Freedom: {df_diff:.0f}")
print(f"Resulting P-value: {p_value}")