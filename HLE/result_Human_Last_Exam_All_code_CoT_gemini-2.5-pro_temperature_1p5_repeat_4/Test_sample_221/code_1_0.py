import pandas as pd
import numpy as np
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm

# Set a random seed for reproducibility, ensuring the result is consistent.
np.random.seed(123)

n = 100
# The R code `1:2` with n=100 effectively generates alternating 1s and 2s.
water_availability = np.tile([1, 2], n // 2)

# Simulate the data based on the provided R script's logic.
df = pd.DataFrame({
    'water_availability': water_availability,
    'plant_intra_div': np.random.choice([1, 2], size=n),
    'plant_inter_div': np.random.choice(range(1, 11), size=n)
})

# Generate the dependent variables according to the specified relationships.
# insect_richness is a function of the initial predictors.
df['insect_richness'] = (df['water_availability'] * 0.01 +
                         df['plant_intra_div'] * 0.5 +
                         df['plant_inter_div'] * 1.2 +
                         np.random.normal(size=n))

# herbivory is a function of insect_richness (with a large coefficient) and other predictors.
df['herbivory'] = (df['insect_richness'] * 3.14 +
                   df['water_availability'] * 0.5 +
                   df['plant_intra_div'] * 0.1 +
                   df['plant_inter_div'] * 0.2 +
                   np.random.normal(size=n))

# Define the formulas for the two models, analogous to the lavaan models.
# Model 1 is the full model, reflecting the data generation process.
model_1_formula = 'herbivory ~ insect_richness + water_availability + plant_intra_div + plant_inter_div'
# Model 2 is the restricted model, which omits the crucial insect_richness term.
model_2_formula = 'herbivory ~ water_availability + plant_intra_div + plant_inter_div'

# Fit the two linear models.
fit_1 = ols(model_1_formula, data=df).fit() # Full model
fit_2 = ols(model_2_formula, data=df).fit() # Restricted model

# Compare the two nested models using an ANOVA.
# This test is analogous to the chi-square difference test performed by lavaan's anova().
# The null hypothesis is that the simpler model (fit_2) is sufficient.
anova_results = anova_lm(fit_2, fit_1)
p_value = anova_results['PR(>F)'][1]

# To fulfill the request to "output each number in the final equation",
# we print the ANOVA table, which contains the components used to calculate
# the final F-statistic and p-value.
print("--- ANOVA Comparison of Nested Models ---")
print(anova_results)

# Below we explicitly show the numbers that go into the F-statistic calculation.
rss_diff = anova_results['sum_sq'][1]
df_diff = anova_results['df'][1]
rss_full_model = anova_results['ss_resid'][1]
df_resid_full_model = anova_results['df_resid'][1]
F_statistic = anova_results['F'][1]

print("\n--- F-Statistic Calculation Breakdown ---")
print(f"Model 1 (Full) Residual Sum of Squares: {rss_full_model:.4f}")
print(f"Model 2 (Restricted) Residual Sum of Squares: {anova_results['ss_resid'][0]:.4f}")
print(f"Difference in Sum of Squares (numerator term): {rss_diff:.4f}")
print(f"Difference in Degrees of Freedom (numerator term): {df_diff:.0f}")
print(f"F-statistic = ({rss_diff:.4f} / {df_diff:.0f}) / ({rss_full_model:.4f} / {df_resid_full_model:.0f}) = {F_statistic:.4f}")
print(f"The resulting p-value for this F-statistic is: {p_value}")