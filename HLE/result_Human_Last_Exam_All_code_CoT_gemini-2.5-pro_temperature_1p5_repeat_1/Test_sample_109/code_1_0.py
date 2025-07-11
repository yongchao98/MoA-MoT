import numpy as np
import statsmodels.api as sm
from statsmodels.sandbox.regression.gmm import IV2SLS

# --- Step 1: Simulate Data ---
# We will simulate data where we have one known valid instrument (Z1)
# and one invalid instrument (Z2), fulfilling the prompt's condition.

np.random.seed(42)  # for reproducibility
n_obs = 1000  # number of observations

# Unobserved confounder U that causes endogeneity (affects both D and Y)
u = np.random.normal(0, 1, n_obs)

# Z1: A valid instrument. It's correlated with D but NOT with Y's error term.
z1 = np.random.normal(0, 1, n_obs)

# Z2: An INvalid instrument. It's correlated with D, but it also has a
# direct effect on Y, violating the exclusion restriction.
z2 = np.random.normal(0, 1, n_obs)

# D: The binary treatment. Its probability depends on both instruments and the confounder.
# We use a latent variable d_star to generate the binary D.
d_star = 1.2 * z1 + 1.2 * z2 + 1 * u + np.random.normal(0, 1, n_obs)
d = (d_star > 0).astype(int)

# Y: The outcome. It depends on the treatment D, the confounder U, and critically,
# it also depends directly on the invalid instrument Z2. The true causal effect of D is 2.
# The coefficient on z2 (0.7) makes it an invalid instrument.
y = 2.0 * d + 0.7 * z2 + 1 * u + np.random.normal(0, 1, n_obs)


# --- Step 2: Perform IV Regression and the Hansen J-test ---
# The model is "overidentified" because we have 2 instruments (Z1, Z2)
# for 1 endogenous variable (D). The J-test checks if the instruments are valid.

# Prepare variables for statsmodels
Y = y
# Add a constant to the model, and specify D as the endogenous variable
X = sm.add_constant(d)
# Add a constant, and specify Z1 and Z2 as the instruments
Z = sm.add_constant(np.column_stack([z1, z2]))

# Fit the Two-Stage Least Squares (2SLS) model
iv_model = IV2SLS(endog=Y, exog=X, instrument=Z).fit()

# The statsmodels result for IV2SLS includes the Sargan-Hansen test of overidentifying restrictions.
# We extract the J-statistic and its corresponding p-value.
sargan_stat = iv_model.j
sargan_p_value = iv_model.j_p


# --- Step 3: Print and Explain the Result ---
print("--- Testing Instrument Validity with Hansen J-Test ---")
print("Scenario: We assume instrument Z1 is valid and want to test the validity of Z2.")
print(f"Model is overidentified: 2 instruments (Z1, Z2) for 1 endogenous variable (D).\n")

print("Null Hypothesis (H0): All instruments (Z1 and Z2) are valid (exogenous).")
print("Alternative Hypothesis (H1): At least one instrument is invalid.\n")

print(f"Hansen J-test statistic: {sargan_stat:.4f}")
print(f"P-value: {sargan_p_value:.4f}\n")

alpha = 0.05
if sargan_p_value < alpha:
    print(f"Result: With a p-value of {sargan_p_value:.4f}, which is less than {alpha}, we REJECT the null hypothesis.")
    print("Conclusion: The test provides strong evidence that not all instruments are valid.")
    print("Since we operated under the assumption that Z1 is valid, we conclude that Z2 is an invalid instrument.")
else:
    print(f"Result: With a p-value of {sargan_p_value:.4f}, we FAIL to reject the null hypothesis.")
    print("Conclusion: The test does not provide evidence to suggest any of the instruments are invalid.")

<<<B>>>