import numpy as np
import pandas as pd
import statsmodels.api as sm
from linearmodels.iv import IV2SLS

# 1. Simulate data
np.random.seed(42)
n = 2000  # Sample size

# Create correlated errors for endogeneity
cov = [[1.0, 0.7], [0.7, 1.0]]  # Covariance matrix for [u, v]
errors = np.random.multivariate_normal([0, 0], cov, n)
u = errors[:, 0]  # Error in outcome equation
v = errors[:, 1]  # Error in treatment equation

# Generate instruments
Z1 = np.random.randn(n)  # Z1 is a valid instrument (exogenous)
# Z2 is an invalid instrument because it's correlated with the outcome's error term 'u'
Z2 = 0.4 * Z1 + 0.6 * u + np.random.randn(n)

# Create the endogenous treatment variable D
# D depends on Z1, Z2, and the error 'v' which is correlated with 'u'
D = 0.5 + 1.2 * Z1 + 0.8 * Z2 + v

# Create the continuous outcome variable Y
# The true effect of D on Y is 2.0
Y = 1.0 + 2.0 * D + u

# 2. Prepare data for the model
# Add a constant for the intercept
X = sm.add_constant(D, prepend=True)
X.columns = ['const', 'D']
Z = sm.add_constant(np.column_stack((Z1, Z2)), prepend=True)
Z.columns = ['const_instr', 'Z1', 'Z2']

# 3. Fit the IV model and perform the test
# We regress Y on D, using Z1 and Z2 as instruments
# The model is overidentified (2 instruments, 1 endogenous regressor)
model = IV2SLS(dependent=Y, exog=X[['const']], endog=X[['D']], instruments=Z[['Z1', 'Z2']])
results = model.fit(cov_type='robust')

# 4. Print results, focusing on the Hansen J-test
print("IV Regression Results")
print("=====================")
print(results)

# The test of overidentifying restrictions is at the bottom of the summary
# It is labeled "Sargan's J" or similar.
j_stat = results.sargan.stat
p_value = results.sargan.pval

print("\n--- Interpreting the Test for Instrument Validity ---")
print(f"Hansen J-statistic: {j_stat:.4f}")
print(f"P-value: {p_value:.4f}")
print("\nNull Hypothesis (H0): All instruments are valid (exogenous).")
print("Alternative Hypothesis (H1): At least one instrument is invalid.")

if p_value < 0.05:
    print("\nConclusion: We REJECT the null hypothesis at the 5% significance level.")
    print("This suggests that at least one instrument is invalid (not exogenous).")
    print("Given our assumption that Z1 is valid, this test provides evidence that Z2 is invalid.")
else:
    print("\nConclusion: We FAIL to reject the null hypothesis.")
    print("There is not enough evidence to conclude that any of the instruments are invalid.")
