import numpy as np
import pandas as pd
from linearmodels.iv import IV2SLS
import statsmodels.api as sm

def run_validity_test():
    """
    This function demonstrates the use of the Hansen J-test to check
    instrument validity.

    Scenario:
    - We have one known valid instrument (Z1).
    - We want to test if a second instrument (Z2) is also valid.
    - We will construct Z2 to be *invalid* to show that the test can detect this.
    """
    # 1. Simulate data
    np.random.seed(42)
    n_obs = 2000

    # Unobserved confounder that affects both treatment and outcome
    u = np.random.normal(0, 1, n_obs)
    
    # Error term for the outcome equation
    v = np.random.normal(0, 1, n_obs)

    # Z1: A valid instrument (uncorrelated with the outcome's error term, v)
    z1 = np.random.normal(0, 1, n_obs)

    # Z2: An INVALID instrument. It violates the exclusion restriction because it
    # is directly correlated with the outcome's error term, v.
    z2 = 0.6 * z1 + 0.8 * v

    # Treatment (D) is determined by the instruments and the confounder
    # A latent variable model is used for the binary treatment
    d_latent = 0.8 * z1 + 0.7 * z2 + 1.2 * u + np.random.normal(0, 1, n_obs)
    D = (d_latent > np.mean(d_latent)).astype(int)

    # Outcome (Y) is determined by the treatment and the confounder
    # True causal effect of D on Y is 2.5
    Y = 2.0 + 2.5 * D + 1.5 * u + v

    # Create a pandas DataFrame
    data = pd.DataFrame({'Y': Y, 'D': D, 'Z1': z1, 'Z2': z2})
    
    print("--- The Hansen J-Test for Overidentifying Restrictions ---")
    print("\nNull Hypothesis (H0): All instruments are valid (exogenous).")
    print("Assumption: We know Z1 is a valid instrument.")
    print("Test: We check the joint validity of {Z1, Z2}. If we reject H0,")
    print("it must be because Z2 is invalid.\n")

    # 2. Perform 2SLS Regression
    # We use both Z1 and Z2 as instruments for the endogenous variable D.
    # The model is overidentified (2 instruments for 1 endogenous regressor).
    exog = sm.add_constant(np.empty((n_obs, 0))) # Only a constant
    endog = data[['D']]
    instruments = data[['Z1', 'Z2']]

    model = IV2SLS(dependent=data['Y'], exog=exog, endog=endog, instruments=instruments)
    results = model.fit(cov_type='robust')

    # 3. Get the Hansen J-test results
    # The '.sargan' attribute in linearmodels performs this test.
    j_test = results.sargan

    j_statistic = j_test.stat
    p_value = j_test.pval
    
    # 4. Print and interpret the results
    print("Hansen J-Statistic:", round(j_statistic, 4))
    print("P-value:", round(p_value, 4))

    print("\n--- Interpretation ---")
    if p_value < 0.05:
        print(f"The p-value ({p_value:.4f}) is less than 0.05.")
        print("We REJECT the null hypothesis that all instruments are valid.")
        print("Since we assumed Z1 is valid, this test provides strong evidence")
        print("that the additional instrument, Z2, is INVALID.")
    else:
        print(f"The p-value ({p_value:.4f}) is not less than 0.05.")
        print("We FAIL to reject the null hypothesis. There is no statistical evidence")
        print("to suggest that the additional instrument, Z2, is invalid.")

# Run the demonstration
run_validity_test()
