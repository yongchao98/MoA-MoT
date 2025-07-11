import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
from sklearn.preprocessing import scale

# Set a seed for reproducibility
np.random.seed(42)

# 1. Simulate the data as described in the R code
n = 100
# Create the initial dataframe
df = pd.DataFrame({
    'water_availability': np.tile([1, 2], n // 2),
    'plant_intra_div': np.random.randint(1, 3, n),
    'plant_inter_div': np.random.randint(1, 11, n)
})

# 2. Generate the dependent variables based on the formulas
# Generate insect_richness
df['insect_richness'] = (df['water_availability'] * 0.01 + 
                         df['plant_intra_div'] * 0.5 + 
                         df['plant_inter_div'] * 1.2 + 
                         np.random.normal(size=n))

# Generate herbivory
df['herbivory'] = (df['insect_richness'] * 3.14 + 
                   df['water_availability'] * 0.5 + 
                   df['plant_intra_div'] * 0.1 + 
                   df['plant_inter_div'] * 0.2 + 
                   np.random.normal(size=n))

# 3. Scale the data (z-score standardization), same as R's scale()
df_scaled = pd.DataFrame(scale(df), columns=df.columns)

# 4. Define and fit the two models, analogous to model_1 and model_2
# This corresponds to the 'herbivory' equation in the SEMs

# Full model (like model_1) includes insect_richness
model_1_fit = smf.ols(
    'herbivory ~ insect_richness + water_availability + plant_intra_div + plant_inter_div',
    data=df_scaled
).fit()

# Reduced model (like model_2) omits insect_richness
model_2_fit = smf.ols(
    'herbivory ~ water_availability + plant_intra_div + plant_inter_div',
    data=df_scaled
).fit()

# 5. Compare the two nested models using ANOVA
# This is analogous to anova(fit_1, fit_2) in lavaan
anova_results = sm.stats.anova_lm(model_2_fit, model_1_fit)

# The p-value is in the 'Pr(>F)' column for the second row
p_value = anova_results.iloc[1]['Pr(>F)']

# 6. Print the resulting P value
# Because the effect of insect_richness is so strong (true coefficient = 3.14),
# the p-value will be extremely small, effectively zero.
print(f"The P value from the model comparison is: {p_value}")
print("This value is extremely close to zero, supporting the choice that the expected p-value is 0.")
