import numpy as np
import pandas as pd
import statsmodels.api as sm

# 0. Set up simulation parameters
N = 100000
np.random.seed(42)
# The treatment effect is +500 for every individual.
true_treatment_effect = 500

# 1. Generate a population with a pre-treatment variable X (e.g., prior income)
# X is centered at 50,000 with a standard deviation of 10,000.
X = np.random.normal(loc=50000, scale=10000, size=N)

# 2. Define potential outcomes Y(0) and Y(1)
# Y(0) (income without treatment) depends on prior income X.
# Y(1) (income with treatment) is Y(0) plus the treatment effect.
# u is random noise representing other factors affecting income.
u = np.random.normal(loc=0, scale=5000, size=N)
Y0 = 10000 + 1.1 * X + u
Y1 = Y0 + true_treatment_effect

# Create a DataFrame to hold our population data
df = pd.DataFrame({'X': X, 'Y0': Y0, 'Y1': Y1})

print("--- Case 1: D is randomly assigned. Regression: Y ~ D ---")
# D is assigned with 50% probability, independent of X.
df['D1'] = np.random.binomial(1, 0.5, size=N)
df['Y1_realized'] = df['D1'] * df['Y1'] + (1 - df['D1']) * df['Y0']
model1 = sm.OLS.from_formula("Y1_realized ~ D1", data=df).fit()
print(f"The estimated regression is: Y = {model1.params['Intercept']:.2f} + {model1.params['D1']:.2f} * D")
print(f"The coefficient on D is {model1.params['D1']:.2f}, which is positive and close to the true effect ({true_treatment_effect}).\n")


print("--- Case 2: D depends on X. Regression: Y ~ D (Omitted Variable) ---")
# D is assigned based on X: lower X means higher probability of treatment.
# This creates a strong negative selection bias.
prob_d2 = 1 / (1 + np.exp((df['X'] - 50000) / 5000))
df['D2'] = np.random.binomial(1, prob_d2, size=N)
df['Y2_realized'] = df['D2'] * df['Y1'] + (1 - df['D2']) * df['Y0']
model2 = sm.OLS.from_formula("Y2_realized ~ D2", data=df).fit()
print(f"The estimated regression is: Y = {model2.params['Intercept']:.2f} + {model2.params['D2']:.2f} * D")
print(f"The coefficient on D is {model2.params['D2']:.2f}. Here it's negative because of strong selection bias, despite a true positive effect.\n")


print("--- Case 3: D depends on X. Regression: Y ~ D + X (Controlled) ---")
# Using the same data as Case 2, but now controlling for X in the regression.
model3 = sm.OLS.from_formula("Y2_realized ~ D2 + X", data=df).fit()
print(f"The estimated regression is: Y = {model3.params['Intercept']:.2f} + {model3.params['D2']:.2f} * D + {model3.params['X']:.2f} * X")
print(f"The coefficient on D is {model3.params['D2']:.2f}, which is positive after controlling for the confounder X.")
