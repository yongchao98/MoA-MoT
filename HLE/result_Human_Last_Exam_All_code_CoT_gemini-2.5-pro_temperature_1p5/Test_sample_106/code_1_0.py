import pandas as pd
import numpy as np
import statsmodels.api as sm

# Set a seed for reproducibility of the random data
np.random.seed(42)

# --- 1. Define Population and Potential Outcomes ---
# Let's create a population of 100,000 individuals
N = 100000
# X is pre-program income, normally distributed around $50k
X = np.random.normal(loc=50000, scale=15000, size=N)
# Income cannot be negative
X = np.maximum(X, 0)

# Y0 is the potential income *without* the jobs program. It depends on prior income X.
# We add some random noise.
u0 = np.random.normal(loc=0, scale=5000, size=N)
Y0 = 20000 + 0.8 * X + u0

# Y1 is the potential income *with* the jobs program.
# We assume a constant positive treatment effect of $8,000 for everyone.
# The assumption "treatment effect for everyone is positive" is satisfied.
ATE = 8000
Y1 = Y0 + ATE

# Create a master DataFrame
df = pd.DataFrame({'X': X, 'Y0': Y0, 'Y1': Y1})

# --- Case 1: D is randomly assigned ---
print("--- Case 1: D is randomly assigned. Regression: Y ~ D ---")
df['D1'] = np.random.binomial(n=1, p=0.5, size=N)
df['Y1_obs'] = df['D1'] * df['Y1'] + (1 - df['D1']) * df['Y0']
model1 = sm.OLS.from_formula("Y1_obs ~ D1", data=df).fit()
# The summary table contains all the numbers for the final equation
print(model1.summary().tables[1])
print(f"Conclusion for Case 1: The coefficient on D1 is {model1.params['D1']:.2f}, which is positive and close to the true ATE of ${ATE}.")
print("-" * 60)


# --- Setup for Case 2 & 3: D is assigned conditional on X ---
# Let's create a selection mechanism where people with lower income (X) are more
# likely to join the program. This will create confounding.
X_std = (df['X'] - df['X'].mean()) / df['X'].std()
prob_d = 1 / (1 + np.exp(1.5 * X_std)) # Lower X -> higher probability
df['D23'] = np.random.binomial(n=1, p=prob_d, size=N)
df['Y23_obs'] = df['D23'] * df['Y1'] + (1 - df['D23']) * df['Y0']

# --- Case 2: Conditional assignment, but regression omits X ---
print("\n--- Case 2: D assigned by X. Regression: Y ~ D (omitting X) ---")
model2 = sm.OLS.from_formula("Y23_obs ~ D23", data=df).fit()
print(model2.summary().tables[1])
print(f"Conclusion for Case 2: The coefficient on D23 is {model2.params['D23']:.2f}, which is negative.")
print("This demonstrates that the coefficient is NOT necessarily positive due to selection bias.")
print("-" * 60)


# --- Case 3: Conditional assignment, regression includes X ---
print("\n--- Case 3: D assigned by X. Regression: Y ~ D + X ---")
model3 = sm.OLS.from_formula("Y23_obs ~ D23 + X", data=df).fit()
print(model3.summary().tables[1])
print(f"Conclusion for Case 3: The coefficient on D23 is {model3.params['D23']:.2f}, which is positive.")
print("By controlling for the confounder X, we recover an estimate close to the true ATE of ${ATE}.")
print("-" * 60)
