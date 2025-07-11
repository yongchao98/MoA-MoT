import numpy as np
import pandas as pd
from linearmodels.iv import IV2SLS
import statsmodels.api as sm

# Step 1: Simulate data according to the user's problem
# We will create a scenario with one valid instrument (z1) and one invalid instrument (z2).
np.random.seed(42)
n_obs = 1000

# z1 is a valid instrument: it will be correlated with D but not directly with Y.
z1 = np.random.normal(0, 1, n_obs)

# z2 is an invalid instrument: it will be correlated with D and also directly with Y,
# violating the exclusion restriction.
z2 = np.random.normal(0, 1, n_obs)

# u is an unobserved confounder that affects both D and Y, making D endogenous.
u = np.random.normal(0, 1, n_obs)

# D is the binary treatment. Its value depends on the instruments and the confounder.
# We use a latent variable model for the binary treatment.
d_latent = 1.5 * z1 + 1.5 * z2 + 1.0 * u + np.random.normal(0, 1, n_obs)
D = (d_latent > 0).astype(int)

# Y is the continuous outcome.
# The true effect of D on Y is 2.0.
# The confounder u affects Y.
# The invalid instrument z2 also directly affects Y (the violation).
Y = 2.0 * D + 1.0 * u + 3.0 * z2 + np.random.normal(0, 1, n_obs)

# Create a pandas DataFrame
data = pd.DataFrame({'Y': Y, 'D': D, 'z1': z1, 'z2': z2})
data = sm.add_constant(data, prepend=True)

# Step 2: Perform IV regression using both instruments
# The model is overidentified because we have 2 instruments (z1, z2) for 1 endogenous variable (D).
# The formula syntax [D ~ z1 + z2] tells the model that D is endogenous and z1, z2 are instruments.
model = IV2SLS.from_formula('Y ~ 1 + [D ~ z1 + z2]', data)
results = model.fit(cov_type='robust')

# Step 3: Perform the Hansen J-test (test of overidentifying restrictions)
# The null hypothesis (H0) is that all instruments are valid (exogenous).
# The alternative (H1) is that at least one instrument is invalid.
# In our setup, we assume z1 is valid, so we are testing the validity of z2.
sargan_test = results.sargan

print("--- Instrumental Variable (IV) Regression Results ---")
print(results)

print("\n--- Test of Overidentifying Restrictions (Hansen J-test) ---")
print(f"This test checks if the instruments are jointly valid.")
print(f"Null Hypothesis (H0): All instruments are valid (exogenous).")
print(f"Since we constructed z2 to be an invalid instrument, we expect to reject H0.")
print("\nTest Results:")
print(f"J-Statistic: {sargan_test.stat:.4f}")
print(f"P-value: {sargan_test.pval:.4f}")

if sargan_test.pval < 0.05:
    print("\nConclusion: The p-value is less than 0.05, so we reject the null hypothesis.")
    print("This suggests that at least one instrument is invalid.")
    print("Given our assumption that z1 is valid, this test correctly indicates that z2 is invalid.")
else:
    print("\nConclusion: The p-value is greater than 0.05, so we fail to reject the null hypothesis.")
    print("We do not have sufficient evidence to claim that any of the instruments are invalid.")
