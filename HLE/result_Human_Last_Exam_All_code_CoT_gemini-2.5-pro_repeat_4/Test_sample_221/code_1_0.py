import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.formula.api import ols

# Plan:
# 1. Recreate the data simulation process from the R code.
# 2. Fit two nested linear models that correspond to the herbivory equations in model_1 and model_2.
# 3. Perform an ANOVA to compare the two nested models, which is analogous to the chi-square difference test in lavaan for this case.
# 4. The test's p-value will show whether removing the 'insect_richness -> herbivory' path significantly worsens the model fit.
# 5. Print the key results of the statistical test.

# Set a seed for reproducibility
np.random.seed(123)

# 1. Simulate the data
n = 100
df = pd.DataFrame({
    # np.tile mimics R's recycling rule for vectors
    'water_availability': np.tile([1, 2], n // 2),
    'plant_intra_div': np.random.randint(1, 3, n),
    'plant_inter_div': np.random.randint(1, 11, n)
})

# Generate dependent variables according to the formulas
df['insect_richness'] = (df['water_availability'] * 0.01 +
                         df['plant_intra_div'] * 0.5 +
                         df['plant_inter_div'] * 1.2 +
                         np.random.randn(n))

# The key relationship: herbivory is strongly dependent on insect_richness
df['herbivory'] = (df['insect_richness'] * 3.14 +
                   df['water_availability'] * 0.5 +
                   df['plant_intra_div'] * 0.1 +
                   df['plant_inter_div'] * 0.2 +
                   np.random.randn(n))

# Scale the data (z-score standardization)
df_scaled = pd.DataFrame(sm.preprocessing.scale(df), columns=df.columns)

# 2. Fit the two nested models
# This analysis focuses on the regression for 'herbivory', as this is the only difference between model_1 and model_2.

# Full model (like in fit_1) includes the path from insect_richness
model_full = ols('herbivory ~ insect_richness + water_availability + plant_intra_div + plant_inter_div', data=df_scaled).fit()

# Reduced model (like in fit_2) omits the path from insect_richness
model_reduced = ols('herbivory ~ water_availability + plant_intra_div + plant_inter_div', data=df_scaled).fit()

# 3. Perform ANOVA to compare the models
# The anova_lm function uses an F-test to compare nested models. The p-value tells us if the more complex model is a significantly better fit.
anova_results = sm.stats.anova_lm(model_reduced, model_full)

# 4. Extract and print the statistics from the comparison
p_value = anova_results.iloc[1]['Pr(>F)']
f_statistic = anova_results.iloc[1]['F']
df_numerator = anova_results.iloc[1]['df_diff']
df_denominator = model_full.df_resid

# 5. Print the final results in a clear format
print("An F-test is used to compare the reduced model (model_2) against the full model (model_1).")
print("The final test result equation is:")
print(f"F(df_numerator={int(df_numerator)}, df_denominator={int(df_denominator)}) = {f_statistic:.2f}, P-value = {p_value}")
print("\nConclusion:")
print("Because the data was generated with a strong effect of 'insect_richness' on 'herbivory', the model that omits this effect (model_2) fits the data very poorly compared to the full model.")
print("The test result shows a highly significant difference, with a P-value that is effectively zero.")
