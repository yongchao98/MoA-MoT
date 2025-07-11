import numpy as np
import pandas as pd
from linearmodels.iv import IV2SLS

# Step 1: Explain the chosen test and logic
print("The correct choice is the Hansen J-test (or Sargan test).")
print("This test assesses the validity of overidentifying restrictions.\n")
print("Null Hypothesis (H0): All instruments are valid (uncorrelated with the outcome equation's error).")
print("Alternative Hypothesis (H1): At least one instrument is invalid.\n")
print("In your scenario, you assume one instrument is valid. If the J-test rejects H0,")
print("it provides evidence that one or more of the *other* instruments must be invalid.\n")
print("---" * 15)
print("Demonstration: Simulating data and running the test.\n")

# Step 2: Simulate data
# We will create a scenario with one valid instrument (Z1) and one invalid one (Z2).
np.random.seed(42)
n_obs = 1000

# u is the error term in the outcome equation.
u = np.random.normal(0, 1, n_obs)

# Z1 is a valid instrument: it's not correlated with u.
z1 = np.random.normal(0, 1, n_obs)

# Z2 is an INvalid instrument: it is correlated with u, violating the exclusion restriction.
z2 = 0.6 * u + np.random.normal(0, 1, n_obs)

# The treatment D is endogenous because it's determined by u (and the instruments).
# We make it binary for this example.
latent_d = 1 + 0.8 * z1 + 0.8 * z2 + 0.5 * u + np.random.normal(0,1,n_obs)
d = (latent_d > np.median(latent_d)).astype(int)

# The outcome Y is determined by the treatment D and the error u.
# The true effect of D on Y is 2.
Y = 10 + 2 * d + u

# Create a pandas DataFrame
data = pd.DataFrame({'Y': Y, 'D': d, 'Z1': z1, 'Z2': z2})
data['const'] = 1

# Step 3: Perform 2SLS Regression and the Hansen J-Test
# Dependent: Y
# Exogenous: A constant
# Endogenous: D
# Instruments: Z1, Z2
model = IV2SLS(dependent=data.Y,
               exog=data['const'],
               endog=data.D,
               instruments=data[['Z1', 'Z2']])

results = model.fit(cov_type='robust')

# Step 4: Print and interpret the results
print("We run a 2SLS regression with Y as the outcome, D as the endogenous treatment,")
print("and Z1 and Z2 as our instruments.")
print("Since Z2 was created to be invalid, we expect the J-test to fail (i.e., have a low p-value).\n")

print(results)

# The Sargan test is the J-test under homoskedasticity.
# The linearmodels library reports it in the summary.
sargan_stat = results.sargan.stat
sargan_p_value = results.sargan.pval

print("\n---" * 15)
print("Explicit Test Result:")
# Final printout showing each number clearly as requested
print(f"Sargan (J-test) Statistic: {sargan_stat:.4f}")
print(f"P-value of Sargan Test: {sargan_p_value:.4f}")

if sargan_p_value < 0.05:
    print("\nConclusion: The p-value is less than 0.05, so we reject the null hypothesis.")
    print("We conclude that at least one instrument is invalid.")
else:
    print("\nConclusion: The p-value is greater than 0.05, so we fail to reject the null hypothesis.")
    print("We cannot conclude that any of the instruments are invalid.")
