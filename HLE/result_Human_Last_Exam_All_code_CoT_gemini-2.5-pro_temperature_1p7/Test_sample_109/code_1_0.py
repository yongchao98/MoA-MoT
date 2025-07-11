import numpy as np
import pandas as pd
from linearmodels.iv import IV2SLS
import statsmodels.api as sm

def test_instrument_validity():
    """
    This function simulates data and performs a Hansen J-test to check for
    instrument validity.

    Scenario:
    - Y: Outcome
    - D: Endogenous binary treatment
    - Z1: A valid instrument
    - Z2: An invalid instrument (violates exclusion restriction)

    The goal is to test the validity of Z2, assuming Z1 is valid.
    The Hansen J-test is the appropriate tool for this.
    """
    # 1. Set up simulation parameters
    np.random.seed(42)
    n_obs = 2000
    true_beta = 2.0  # True effect of treatment D on outcome Y

    # 2. Generate data
    # Generate the instruments
    Z1 = np.random.randn(n_obs)  # Valid instrument
    Z2 = np.random.randn(n_obs)  # Instrument to be tested

    # Generate correlated errors to create endogeneity
    # e_u is the error in the outcome equation
    # e_d is the error in the treatment equation
    cov_matrix = [[1.0, 0.6], [0.6, 1.0]]
    errors = np.random.multivariate_normal([0, 0], cov_matrix, size=n_obs)
    e_u = errors[:, 0]
    e_d = errors[:, 1]

    # Create the structural error 'u' for the outcome equation.
    # We make Z2 invalid by having it directly influence the error term 'u'.
    # This violates the exclusion restriction (Cov(Z2, u) != 0).
    # Z1 remains valid as Cov(Z1, u) = 0 by construction.
    u = 0.8 * Z2 + e_u

    # Create the treatment variable D.
    # It depends on both Z1 and Z2 (relevance condition).
    # It also depends on e_d, which is correlated with u, making D endogenous.
    d_latent = 1.0 * Z1 + 1.0 * Z2 + e_d
    D = (d_latent > 0).astype(int)

    # Create the outcome variable Y
    Y = 10.0 + true_beta * D + u

    # 3. Create a pandas DataFrame for easier handling
    df = pd.DataFrame({'Y': Y, 'D': D, 'Z1': Z1, 'Z2': Z2})
    # Add a constant for the regression model
    df['const'] = 1

    # 4. Perform IV (2SLS) regression and conduct the J-test
    # We regress Y on D, using Z1 and Z2 as instruments.
    # The model is overidentified (2 instruments, 1 endogenous regressor).
    model = IV2SLS(dependent=df['Y'],
                   exog=df['const'],
                   endog=df['D'],
                   instruments=df[['Z1', 'Z2']])

    results = model.fit(cov_type='robust')

    # 5. Print and interpret the results
    print("--- IV Regression Results ---")
    print(results)
    print("\n--- Test of Overidentifying Restrictions (Hansen J-Test) ---")
    
    # The J-statistic is available in the results object
    j_stat = results.j_stat.stat
    j_pval = results.j_stat.pval

    print(f"Hansen J-statistic: {j_stat:.4f}")
    print(f"P-value of J-statistic: {j_pval:.4f}")

    print("\n--- Interpretation ---")
    print("The null hypothesis (H0) for the J-test is that all instruments are valid (exogenous).")
    
    if j_pval < 0.05:
        print("Result: The p-value is less than 0.05, so we REJECT the null hypothesis.")
        print("Conclusion: There is strong evidence that at least one instrument is invalid.")
        print("Given our assumption that Z1 is valid, we conclude that Z2 is an invalid instrument.")
    else:
        print("Result: The p-value is greater than 0.05, so we FAIL to reject the null hypothesis.")
        print("Conclusion: We do not have enough evidence to claim that any of the instruments are invalid.")

# Run the demonstration
test_instrument_validity()