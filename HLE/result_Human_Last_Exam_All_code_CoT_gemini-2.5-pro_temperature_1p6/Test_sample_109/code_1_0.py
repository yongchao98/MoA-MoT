import numpy as np
import pandas as pd
import statsmodels.api as sm
from linearmodels.iv import IV2SLS

# --- 1. Explanation ---
# The Hansen J-test (or Sargan-Hansen test) is designed to test the validity of "overidentifying restrictions".
# In an instrumental variable (IV) context, this means you have more instruments than you need to identify the causal effect.
# The null hypothesis (H0) of the J-test is that ALL instruments are valid (i.e., uncorrelated with the error term).
# If the p-value is small, we reject H0, concluding that AT LEAST ONE instrument is invalid.
#
# In your case, you assume one instrument is valid. If you run a J-test with your known valid instrument plus
# the instruments you want to test, a rejection of the null implies that one or more of the "other"
# instruments must be invalid.

# --- 2. Simulation ---
# Let's simulate data where Z1 is a valid instrument and Z2 is invalid.

# Set a seed for reproducibility
np.random.seed(42)

# Sample size
n_obs = 1000

# U is an unobserved confounder that affects both the treatment D and the outcome Y.
# This is what makes D endogenous.
U = np.random.normal(0, 1, n_obs)

# Z1 is a VALID instrument. It affects D but is NOT related to U.
z_valid = np.random.normal(0, 2, n_obs)

# Z2 is an INVALID instrument. It affects D but IS ALSO related to U (violates exogeneity).
z_invalid = 0.6 * U + np.random.normal(0, 1, n_obs)

# The binary treatment D is determined by the instruments and the confounder.
d_latent = 1 + 1.5 * z_valid + 1.5 * z_invalid - U + np.random.normal(0, 5, n_obs)
D = (d_latent > 0).astype(int)

# The continuous outcome Y is determined by the treatment D and the confounder U.
# The true causal effect of D on Y is 2.5.
Y = 10 + 2.5 * D + 2 * U + np.random.normal(0, 2, n_obs)

# Create a pandas DataFrame
data = pd.DataFrame({'Y': Y, 'D': D, 'z_valid': z_valid, 'z_invalid': z_invalid})
data = sm.add_constant(data, prepend=False) # Add a constant for the regression

# --- 3. Perform 2SLS and the J-Test ---
# We use both z_valid and z_invalid as instruments for the endogenous variable D.
# The model is overidentified because there are 2 instruments for 1 endogenous variable.

model = IV2SLS(dependent=data.Y,
               exog=data.const,
               endog=data.D,
               instruments=data[['z_valid', 'z_invalid']])

results = model.fit(cov_type='robust')

# --- 4. Print the Results ---
# The summary includes the Sargan statistic, which is the J-test for overidentifying restrictions.
print("--- 2SLS Estimation Results ---")
print(results)
print("\n--- Interpretation of the J-Test ---")
# The Sargan test is the J-test. Its null hypothesis is that the instruments are valid.
j_statistic = results.sargan.stat
p_value = results.sargan.pval
print(f"Sargan (J-test) Statistic: {j_statistic:.4f}")
print(f"P-value: {p_value:.4f}")

if p_value < 0.05:
    print("\nThe p-value is less than 0.05. We reject the null hypothesis.")
    print("This suggests that at least one instrument is invalid.")
    print("Given our assumption that z_valid is valid, this test provides evidence against the validity of z_invalid.")
else:
    print("\nThe p-value is greater than 0.05. We fail to reject the null hypothesis.")
    print("The test does not provide evidence to suggest any of the instruments are invalid.")
