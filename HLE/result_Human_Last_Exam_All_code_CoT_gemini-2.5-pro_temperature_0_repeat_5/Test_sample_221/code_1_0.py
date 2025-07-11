import pandas as pd
import numpy as np
from statsmodels.formula.api import ols
from statsmodels.stats.api import anova_lm

# Set a seed for reproducibility of the random data
np.random.seed(123)

# 1. Data Generation
n = 100
# Replicate R's recycling of '1:2' to length n
water_availability = np.tile([1, 2], n // 2)
# Generate the other predictors
plant_intra_div = np.random.choice([1, 2], n, replace=True)
plant_inter_div = np.random.choice(np.arange(1, 11), n, replace=True)

# Create a DataFrame for the raw (unscaled) data
example_raw = pd.DataFrame({
    'water_availability': water_availability,
    'plant_intra_div': plant_intra_div,
    'plant_inter_div': plant_inter_div
})

# Generate the dependent variables based on the specified formulas
example_raw['insect_richness'] = (example_raw['water_availability'] * 0.01 +
                                  example_raw['plant_intra_div'] * 0.5 +
                                  example_raw['plant_inter_div'] * 1.2 +
                                  np.random.randn(n))

example_raw['herbivory'] = (example_raw['insect_richness'] * 3.14 +
                            example_raw['water_availability'] * 0.5 +
                            example_raw['plant_intra_div'] * 0.1 +
                            example_raw['plant_inter_div'] * 0.2 +
                            np.random.randn(n))

# Scale the entire dataframe, as done in the R script
example_scaled = (example_raw - example_raw.mean()) / example_raw.std()

# 2. Model Fitting
# Fit the full model (model_1)
fit_1 = ols('herbivory ~ insect_richness + water_availability + plant_intra_div + plant_inter_div', data=example_scaled).fit()

# Fit the restricted model (model_2)
fit_2 = ols('herbivory ~ water_availability + plant_intra_div + plant_inter_div', data=example_scaled).fit()

# 3. Model Comparison using ANOVA
# This test is analogous to lavaan's anova() for nested models.
# It tests the null hypothesis that the simpler model (fit_2) is sufficient.
anova_results = anova_lm(fit_2, fit_1)

# The p-value is in the 'PR(>F)' column for the second model comparison row.
p_value = anova_results['PR(>F)'][1]

# The p-value will be extremely small, effectively 0.
# We print the numbers from the ANOVA test result to show the full context.
# The final equation here is the test itself: F(df_diff, df_resid) = F-statistic, p = p-value
print("ANOVA Comparison of Models")
print("--------------------------")
# The difference in degrees of freedom (number of constrained parameters)
df_diff = int(anova_results['df_diff'][1])
# The residual degrees of freedom for the full model
df_resid = int(anova_results['df_resid'][1])
# The F-statistic of the comparison
f_stat = anova_results['F'][1]

print(f"Difference in DF: {df_diff}")
print(f"Residual DF (Full Model): {df_resid}")
print(f"F-statistic: {f_stat:.4f}")
print(f"P-value: {p_value}")
