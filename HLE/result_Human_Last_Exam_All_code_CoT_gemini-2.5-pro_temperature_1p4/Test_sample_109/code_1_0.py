import numpy as np
import pandas as pd
from linearmodels.iv import IV2SLS

def perform_hansen_j_test():
    """
    This function simulates data for an IV setup and performs the Hansen J-test
    to check for the validity of overidentifying restrictions.
    """
    # 1. Set up the simulation parameters
    n = 1000  # Number of observations
    np.random.seed(42) # for reproducibility
    beta1 = 0.5  # True causal effect of D on Y

    # 2. Generate the data
    # u is the unobserved confounder that affects both treatment and outcome
    u = np.random.normal(0, 1, n)

    # Z1 is a valid instrument: correlated with D but not with u
    z1 = np.random.normal(0, 1, n)

    # Z2 is an invalid instrument: correlated with D AND with u
    # This violates the exogeneity assumption.
    z2 = 0.6 * u + np.random.normal(0, 1, n)

    # Create the endogenous binary treatment D
    # A latent variable d_star determines D. d_star is affected by z1, z2, and u.
    d_star = 0.8 * z1 + 0.8 * z2 + u + np.random.normal(0, 1, n)
    D = (d_star > 0).astype(int)

    # Create the outcome Y
    # Y is affected by D and the confounder u.
    Y = 1.0 + beta1 * D + 1.5 * u + np.random.normal(0, 1, n)

    # 3. Prepare data for the model
    data = pd.DataFrame({'Y': Y, 'D': D, 'Z1': z1, 'Z2': z2})
    data['const'] = 1
    
    # Define dependent, exogenous, endogenous, and instrument variables
    dependent = data['Y']
    exog = data['const']
    endog = data['D']
    instruments = data[['Z1', 'Z2']]

    # 4. Perform the IV regression and Hansen J-test
    # We have more instruments (2) than endogenous variables (1), so the model is overidentified.
    # The IV2SLS function from linearmodels automatically computes the test.
    print("Running Instrumental Variable regression with two instruments (Z1, Z2) for one endogenous regressor (D).")
    print("Assumption: We know Z1 is valid and want to test if Z2 is also valid.")
    print("We use the Hansen J-test (also known as the Sargan test of overidentifying restrictions).\n")

    iv_model = IV2SLS(dependent, exog, endog, instruments)
    iv_results = iv_model.fit(cov_type='robust')

    # The test is available in the results object
    sargan_test = iv_results.sargan

    # 5. Print and interpret the results
    print("--- Hansen J-Test Results ---")
    print(f"Null Hypothesis: All instruments are valid (exogenous).")
    print(f"Test Statistic: {sargan_test.stat:.4f}")
    print(f"P-value: {sargan_test.pval:.4f}")
    print(f"Degrees of Freedom: {sargan_test.df}")
    
    print("\n--- Interpretation ---")
    if sargan_test.pval < 0.05:
        print("The p-value is less than 0.05, so we reject the null hypothesis.")
        print("This suggests that at least one instrument is invalid.")
        print("Given our prior knowledge that Z1 is valid, this test provides strong evidence that Z2 is invalid.")
    else:
        print("The p-value is not less than 0.05, so we fail to reject the null hypothesis.")
        print("The test does not provide evidence to suggest any of the instruments are invalid.")

perform_hansen_j_test()
<<<B>>>