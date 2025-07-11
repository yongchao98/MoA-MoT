import numpy as np
import statsmodels.api as sm
from scipy import stats

# Set a seed for reproducibility, similar to R's default behavior
np.random.seed(123)

n = 100
# Create the dataframe components
water_availability = np.tile(np.array([1, 2]), n // 2)
plant_intra_div = np.random.randint(1, 3, n)
plant_inter_div = np.random.randint(1, 11, n)

# Create the dependent variables based on the specified relationships
insect_richness = (water_availability * 0.01 + 
                   plant_intra_div * 0.5 + 
                   plant_inter_div * 1.2 + 
                   np.random.normal(size=n))

herbivory = (insect_richness * 3.14 + 
             water_availability * 0.5 + 
             plant_intra_div * 0.1 + 
             plant_inter_div * 0.2 + 
             np.random.normal(size=n))

# Combine into a dictionary and scale the data
data = {
    'herbivory': herbivory,
    'insect_richness': insect_richness,
    'water_availability': water_availability,
    'plant_intra.div': plant_intra_div,
    'plant_inter.div': plant_inter_div
}
scaled_data = {key: stats.zscore(value) for key, value in data.items()}

# Define predictors for the models
# This is analogous to model_2 (the restricted model)
X_model_2 = np.column_stack([
    scaled_data['water_availability'],
    scaled_data['plant_intra.div'],
    scaled_data['plant_inter.div']
])
X_model_2 = sm.add_constant(X_model_2) # Add intercept

# This is analogous to model_1 (the full model)
X_model_1 = np.column_stack([
    scaled_data['insect_richness'],
    scaled_data['water_availability'],
    scaled_data['plant_intra.div'],
    scaled_data['plant_inter.div']
])
X_model_1 = sm.add_constant(X_model_1) # Add intercept

y = scaled_data['herbivory']

# Fit the two models
fit_1 = sm.OLS(y, X_model_1).fit()
fit_2 = sm.OLS(y, X_model_2).fit()

# Perform an ANOVA to compare the two nested models
# This is the equivalent test to lavaan's anova for these models
anova_results = sm.stats.anova_lm(fit_2, fit_1)
p_value = anova_results.iloc[1]['Pr(>F)']

# The F-test here tests the null hypothesis that the coefficient for insect_richness is 0.
# Because the data was generated with a very large coefficient (3.14),
# we expect a very large F-statistic and a very small p-value.
print(f"The ANOVA comparison between the full model (model_1) and the reduced model (model_2) tests if the added predictor ('insect_richness') is significant.")
print(f"The null hypothesis is that the coefficient for 'insect_richness' is zero.")
print(f"From the data generation (herbivory = insect_richness * 3.14 + ...), we know the true coefficient is 3.14.")
print(f"Therefore, we expect to strongly reject the null hypothesis, resulting in a very small p-value.")
print(f"\nThe calculated P value is: {p_value}")
print(f"As expected, the p-value is extremely small, effectively zero.")
