import numpy as np
import pandas as pd
from linearmodels.iv import IV2SLS
import statsmodels.api as sm

def run_hansen_j_test_simulation():
    """
    This function demonstrates the use of the Hansen J-test to check
    the validity of surplus instruments.
    """
    # 1. Explanation of the test
    print("--------------------------------------------------------------------")
    print("Hansen J-Test for Overidentifying Restrictions")
    print("--------------------------------------------------------------------")
    print("The Hansen J-test is used when you have more instruments than endogenous variables.")
    print("The null hypothesis (H0) is that all instruments are valid (exogenous).\n")
    print("If we assume one instrument is valid, rejecting H0 implies that the additional instruments are invalid.")
    print("We will simulate data where Z1 is a valid instrument and Z2 is invalid.\n")

    # 2. Simulate Data
    np.random.seed(42)  # for reproducibility
    n_obs = 2000

    # Z1 is a valid instrument (uncorrelated with the outcome error U)
    z1 = np.random.randn(n_obs)
    
    # U is the error term in the outcome equation
    u = np.random.randn(n_obs)

    # Z2 is an invalid instrument because it is correlated with U
    # This violates the exogeneity/exclusion restriction.
    z2 = 0.6 * u + np.random.randn(n_obs)

    # D is the binary endogenous treatment. Its value depends on Z1 and Z2.
    # We create a latent variable D_star and then threshold it.
    d_latent = 1.2 * z1 + 1.0 * z2 + np.random.randn(n_obs)
    d = (d_latent > 0).astype(int)

    # Y is the continuous outcome. It depends on D and the error U.
    # True effect of D on Y is 2.5
    y = 1.5 + 2.5 * d + u

    # Prepare data for the model
    data = pd.DataFrame({'Y': y, 'D': d, 'Z1': z1, 'Z2': z2})
    data['const'] = 1

    # 3. Run 2SLS Regression and Hansen J-Test
    # We regress Y on D, using both Z1 and Z2 as instruments.
    # The model is overidentified (2 instruments, 1 endogenous variable).
    model = IV2SLS(
        dependent=data['Y'],
        exog=data['const'],
        endog=data['D'],
        instruments=data[['Z1', 'Z2']]
    )
    results = model.fit(cov_type='robust')

    # 4. Print and Interpret Results
    print("Running 2SLS with both Z1 (valid) and Z2 (invalid) as instruments...")
    print(results)

    # Extract the J-statistic and its p-value
    j_statistic = results.j_stat.stat
    j_p_value = results.j_stat.pval

    print("\n--------------------------------------------------------------------")
    print("Interpretation of the J-Test Result")
    print("--------------------------------------------------------------------")
    print(f"The Hansen J-statistic is: {j_statistic:.4f}")
    print(f"The p-value of the J-statistic is: {j_p_value:.4f}")

    if j_p_value < 0.05:
        print("\nResult: The p-value is less than 0.05, so we REJECT the null hypothesis.")
        print("Conclusion: There is strong evidence that at least one instrument is invalid.")
        print("Since we assumed Z1 was valid, we conclude that Z2 is an invalid instrument.")
    else:
        print("\nResult: The p-value is greater than 0.05, so we FAIL to reject the null hypothesis.")
        print("Conclusion: We do not have enough evidence to claim that any of the instruments are invalid.")

    # Output the final equation as requested
    b_const = results.params['const']
    b_d = results.params['D']
    print("\n--------------------------------------------------------------------")
    print("Final Estimated Equation")
    print("--------------------------------------------------------------------")
    print(f"Y = {b_const:.4f} + {b_d:.4f} * D")
    print("Note: The coefficient for D is likely biased because we used an invalid instrument.")
    print("--------------------------------------------------------------------")

# Execute the function
run_hansen_j_test_simulation()