import numpy as np
import pandas as pd
from linearmodels.iv import IVGMM

# Step 1: Simulate data according to the problem description
np.random.seed(42)
n_obs = 2000

# Z1 and Z2 are our valid instruments.
Z1 = np.random.normal(size=n_obs)
Z2 = np.random.normal(size=n_obs)

# Z3 is our potentially invalid instrument.
Z3 = np.random.normal(size=n_obs)

# U is an unobserved confounder that makes treatment D endogenous.
U = np.random.normal(size=n_obs)

# The binary treatment D is determined by the instruments and the confounder.
# All instruments are relevant (correlated with D).
d_star = 1.5 * Z1 + 1.5 * Z2 + 1.5 * Z3 + 2 * U + np.random.normal(size=n_obs)
D = (d_star > 0).astype(int)

# The outcome Y is determined by the treatment D and the confounder U.
# CRUCIALLY, Y is also directly affected by Z3, making Z3 an invalid instrument
# because it violates the exclusion restriction.
# The true causal effect of D on Y is 2.0.
Y = 1.0 + 2.0 * D + 3.0 * U + 2.5 * Z3 + np.random.normal(size=n_obs)

# Create a DataFrame
data = pd.DataFrame({'Y': Y, 'D': D, 'Z1': Z1, 'Z2': Z2, 'Z3': Z3})
data['const'] = 1

# Step 2: Perform the Hansen J-test
# We use all three instruments (Z1, Z2, Z3) to instrument for the single
# endogenous variable (D). This is an overidentified model (3 instruments > 1 endogenous var).
# The null hypothesis of the J-test is that all instruments are valid (exogenous).
# Since we constructed Z3 to be invalid, we expect the test to reject the null.

# Define model variables
dependent = data['Y']
exog = data['const']      # Exogenous regressors in the main equation
endog = data['D']         # Endogenous regressor
instruments = data[['Z1', 'Z2', 'Z3']] # Instruments

# Fit the IV-GMM model
model = IVGMM(dependent, exog, endog, instruments).fit(cov_type='robust')

# Extract the J-test results
j_statistic = model.j_stat.stat
p_value = model.j_stat.pval
degrees_of_freedom = model.j_stat.df

# Step 3: Print the results and conclusion
print("Hansen J-Test for Overidentifying Restrictions")
print("----------------------------------------------")
print("This test evaluates the null hypothesis that all instruments are valid (exogenous).")
print("Our model has 1 endogenous variable (D) and 3 instruments (Z1, Z2, Z3).")
print(f"This gives us {degrees_of_freedom} degrees of freedom (num instruments - num endogenous vars).")
print("\nIn our simulated data, instrument Z3 was intentionally made invalid.")
print("Therefore, we expect the test to reject the null hypothesis, indicated by a small p-value.")
print("\nTest Results:")
print(f"  Hansen J-statistic: {j_statistic:.4f}")
print(f"  P-value: {p_value:.4f}")

if p_value < 0.05:
    print("\nConclusion: The p-value is less than 0.05, so we reject the null hypothesis.")
    print("This correctly indicates that at least one of our instruments is invalid.")
else:
    print("\nConclusion: The p-value is greater than 0.05, so we fail to reject the null hypothesis.")
    print("This would incorrectly suggest that our instruments are valid.")
