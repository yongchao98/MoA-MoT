import numpy as np
import pandas as pd
from linearmodels.iv import IV2SLS

# Step 1: Explain the methodology and set up the simulation
print("The correct test is the Hansen J-test (or test of overidentifying restrictions).")
print("This test assesses whether the instruments are jointly valid under the assumption that there are more instruments than needed.")
print("The null hypothesis is that all instruments are valid. If we reject the null, it implies at least one instrument is invalid.")
print("\n--- Simulating Data ---")
print("We create two instruments: Z1 (valid) and Z2 (invalid).")
print("Z2 is invalid because we will make it directly correlated with the error term in the outcome equation for Y.")

# for reproducibility
np.random.seed(42)

# number of observations
n = 1000

# Step 2: Generate the data according to the scenario
# Z1 is a valid instrument
z1 = np.random.normal(0, 1, n)

# Z2 will be our invalid instrument
z2 = np.random.normal(0, 1, n)

# The error term for the outcome equation (Y)
# We make it correlated with Z2, which violates the exogeneity/exclusion condition for Z2.
u = 0.8 * z2 + np.random.normal(0, 1, n)

# The error term for the treatment equation (D)
v = np.random.normal(0, 1, n)

# The binary treatment D is determined by both instruments (relevance condition is met for both)
# The latent variable d_star determines the binary D
d_star = 1.0 * z1 + 1.2 * z2 + v
D = (d_star > 0).astype(int)

# The continuous outcome Y is determined by the treatment D and the error u
# Since u is correlated with Z2, Z2 has an invalid, direct effect on Y
beta = 2.5
Y = 1.5 + beta * D + u

# Create a DataFrame
df = pd.DataFrame({'Y': Y, 'D': D, 'Z1': z1, 'Z2': z2})
df['const'] = 1

# Step 3: Run the IV regression using both instruments
# The model is "overidentified" because we have 2 instruments for 1 endogenous variable.
print("\n--- Performing IV Regression with both Z1 and Z2 as Instruments ---")
iv_model = IV2SLS.from_formula('Y ~ 1 + [D ~ Z1 + Z2]', data=df).fit(cov_type='robust')

# Step 4: Print the results and interpret the Hansen J-test
print("\n--- Regression Results ---")
print(iv_model)

print("\n--- Interpretation of the Test ---")
# The test is reported as "Sargan's J" in the linearmodels output
j_statistic = iv_model.sargan.stat
p_value = iv_model.sargan.pval

print(f"The test for overidentifying restrictions (Hansen/Sargan J-test) gives a J-statistic of: {j_statistic:.4f}")
print(f"The p-value for this test is: {p_value:.4f}")

if p_value < 0.05:
    print("\nSince the p-value is less than 0.05, we REJECT the null hypothesis.")
    print("We conclude that at least one of the instruments (Z1 or Z2) is invalid.")
    print("In our simulation, we know this is correct because we deliberately made Z2 invalid.")
else:
    print("\nSince the p-value is greater than 0.05, we FAIL to reject the null hypothesis.")
