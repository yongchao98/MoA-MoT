import numpy as np
import pandas as pd
from linearmodels.iv import IV2SLS
import statsmodels.api as sm

# Step 1: Plan
# We will simulate data to demonstrate the Hansen J-test for instrument validity.
# 1. Define an unobserved confounder 'U' that affects both the treatment 'D' and the outcome 'Y'.
# 2. Create a valid instrument 'Z1' that is correlated with 'D' but not 'U'.
# 3. Create an invalid instrument 'Z2' that is correlated with 'D' AND 'U' (violating the exclusion restriction).
# 4. Generate the treatment 'D' and outcome 'Y' based on these variables.
# 5. Run a 2-Stage Least Squares (2SLS) regression using both Z1 and Z2 as instruments.
# 6. Examine the Hansen J-statistic. Since one instrument is invalid, we expect the test to reject the null hypothesis that all instruments are valid.

# Step 2: Simulate Data
np.random.seed(42)
n = 1000  # Sample size

# Unobserved confounder
U = np.random.normal(0, 1, n)

# Z1 is a valid instrument: correlated with D (through its own variable), not with U
Z1 = np.random.normal(0, 1, n)

# Z2 is an INvalid instrument: we make it correlated with U by construction
Z2 = 0.7 * U + np.random.normal(0, 1, n)

# First stage: D is determined by Z1, Z2, and the confounder U
# We make D binary for this example
d_latent = 1.2 * Z1 + 0.9 * Z2 + 1.5 * U + np.random.normal(0, 1, n)
D = (d_latent > np.mean(d_latent)).astype(int)

# Second stage (Structural Equation): Y is determined by D and the confounder U
true_beta = 5.0
Y = 10.0 + true_beta * D + 2.5 * U + np.random.normal(0, 2, n)

# Create a pandas DataFrame
data = pd.DataFrame({'Y': Y, 'D': D, 'Z1': Z1, 'Z2': Z2})
data['Intercept'] = 1

# Step 3: Perform 2SLS and the J-test
# We have 1 endogenous variable (D) and 2 instruments (Z1, Z2), so the model is overidentified.
# The null hypothesis of the J-test is that BOTH Z1 and Z2 are valid instruments.
# We are assuming Z1 is valid and want to test if Z2 is also valid.
model = IV2SLS.from_formula('Y ~ 1 + [D ~ Z1 + Z2]', data).fit()

# Step 4: Print the results and interpret
print("Demonstration of the Hansen J-Test")
print("=" * 35)
print(f"We simulated data where 'Z1' is a valid instrument and 'Z2' is invalid.\n")
print("We run a 2SLS regression of Y on D, using both Z1 and Z2 as instruments.")
print("The model is overidentified, allowing us to perform the J-test.\n")
print("Null Hypothesis of the J-test: All instruments (Z1, Z2) are valid (exogenous).")
print("Our assumption: We know 'Z1' is valid.\n")
print("Therefore, if we reject the null, we conclude that 'Z2' must be the invalid instrument.\n")

# Extracting and printing the J-statistic and its p-value
j_stat = model.sargan.stat
j_pval = model.sargan.pval

print("-" * 35)
print("Test Results:")
print(f"Hansen J-statistic: {j_stat:.4f}")
print(f"P-value of J-test: {j_pval:.4f}")
print("-" * 35)

print("\nConclusion:")
if j_pval < 0.05:
    print(f"The p-value ({j_pval:.4f}) is less than 0.05.")
    print("We REJECT the null hypothesis. This provides strong evidence that not all instruments are valid.")
    print("Since we assumed Z1 was valid, we conclude that Z2 is an invalid instrument.")
else:
    print(f"The p-value ({j_pval:.4f}) is greater than 0.05.")
    print("We FAIL TO REJECT the null hypothesis. We do not have evidence to claim Z2 is invalid.")

# Optional: printing the full model summary for context
print("\nFull Model Summary from 'linearmodels':")
print(model)
