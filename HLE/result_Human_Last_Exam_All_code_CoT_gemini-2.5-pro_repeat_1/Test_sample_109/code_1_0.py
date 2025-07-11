import numpy as np
import pandas as pd
from linearmodels.iv import IV2SLS

# Plan:
# 1. Simulate data with one known valid instrument (Z1) and one invalid instrument (Z2).
#    An instrument is invalid if it violates the exclusion restriction, meaning it is
#    directly correlated with the error term of the outcome (Y) equation.
# 2. Run a 2-Stage Least Squares (2SLS) regression of Y on D, using both Z1 and Z2 as instruments.
# 3. Perform the Hansen J-test (also known as the Sargan test or test of overidentifying restrictions).
# 4. Interpret the result. A rejection (small p-value) indicates that at least one
#    instrument is invalid. Given our assumption that Z1 is valid, this would implicate Z2.

# Step 1: Simulate data
np.random.seed(42)  # for reproducibility
n_obs = 2000  # number of observations

# Z1 is our known valid instrument
z1 = np.random.normal(0, 1, n_obs)

# Z2 is the instrument we want to test. We will construct it to be *invalid*.
z2 = np.random.normal(0, 1, n_obs)

# 'u' is the error term in the final outcome equation for Y.
# To make Z2 invalid, we create a correlation between Z2 and u.
# The true relationship is u = 0.7 * Z2 + noise.
# This violates the exclusion restriction, as Z2 now directly affects Y via u.
noise_u = np.random.normal(0, 1, n_obs)
u = 0.7 * z2 + noise_u

# 'v' is the error term for the first-stage equation (predicting D).
noise_v = np.random.normal(0, 1, n_obs)
# Create the endogenous binary treatment D. Both Z1 and Z2 are relevant instruments (correlated with D).
d_latent = 0.5 + 1.2 * z1 + 1.0 * z2 + noise_v
d = (d_latent > 0).astype(int)

# Create the continuous outcome Y. The true causal effect of D on Y is 2.0.
y = 10.0 + 2.0 * d + u

# Assemble the data into a pandas DataFrame
data = pd.DataFrame({'Y': y, 'D': d, 'Z1': z1, 'Z2': z2})
data['const'] = 1 # Add a constant for the intercept

# Step 2: Define and run the IV model
# We are modeling Y as a function of D, with Z1 and Z2 as instruments.
dependent = data['Y']
exog = data['const']      # Exogenous regressors in the final equation (intercept)
endog = data['D']         # Endogenous regressor
instruments = data[['Z1', 'Z2']]  # All instruments to be used

# Fit the IV2SLS model
iv_model = IV2SLS(dependent, exog, endog, instruments)
iv_results = iv_model.fit(cov_type='robust')

# Step 3: Perform and print the Hansen J-test results
sargan_test = iv_results.sargan
sargan_stat = sargan_test.stat
sargan_pvalue = sargan_test.pval
beta_0_hat = iv_results.params['const']
beta_1_hat = iv_results.params['D']

# Step 4: Print all results and interpretation
print("="*60)
print("Instrumental Variable (2SLS) Regression Results")
print("="*60)
print(iv_results)

print("\n" + "="*60)
print("Estimated Final Equation")
print("="*60)
print(f"The estimated equation for the outcome Y is:")
print(f"Y = {beta_0_hat:.4f} + {beta_1_hat:.4f} * D")

print("\n" + "="*60)
print("Test of Overidentifying Restrictions (Hansen J-test)")
print("="*60)
print("Null Hypothesis (H0): All instruments are valid (exogenous).")
print("In this scenario, since we assume Z1 is valid, a rejection of H0")
print("implies that the additional instrument, Z2, must be invalid.")
print("-" * 30)
print(f"Hansen J-statistic: {sargan_stat:.4f}")
print(f"P-value: {sargan_pvalue:.4f}")
print("-" * 30)

if sargan_pvalue < 0.05:
    print("\nConclusion: We REJECT the null hypothesis at the 5% significance level.")
    print("The test provides strong evidence that at least one instrument is invalid.")
    print("Given our premise, we conclude that Z2 is not a valid instrument.")
else:
    print("\nConclusion: We FAIL to reject the null hypothesis at the 5% significance level.")
    print("The test does not provide enough evidence to say that any instrument is invalid.")
