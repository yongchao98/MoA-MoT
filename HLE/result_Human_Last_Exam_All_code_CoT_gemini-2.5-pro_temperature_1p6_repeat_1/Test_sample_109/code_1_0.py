import numpy as np
import pandas as pd
from linearmodels.iv import IV2SLS

def perform_instrument_validity_test():
    """
    This function demonstrates how to use the Hansen J-test to test the
    validity of a subset of instruments, assuming one is valid.
    """
    # Step 1: Simulate data
    # We create a scenario with one valid instrument (Z1) and one invalid instrument (Z2).
    # An instrument is invalid if it violates the exogeneity assumption, i.e.,
    # it is correlated with the error term 'u' in the outcome equation.
    np.random.seed(42)
    n_obs = 1000

    # Unobserved confounder (becomes the error term in the model)
    u = np.random.normal(0, 1, n_obs)

    # Z1: The known VALID instrument (uncorrelated with u)
    z1 = np.random.normal(0, 1, n_obs)

    # Z2: The suspect INVALID instrument (we make it correlated with u)
    z2 = 0.6 * u + np.random.normal(0, np.sqrt(1 - 0.6**2), n_obs)

    # D: The binary treatment. It's endogenous because it's influenced by u.
    # It also depends on the instruments to satisfy the relevance condition.
    d_latent = 1.0 * z1 + 1.0 * z2 + 0.7 * u + np.random.normal(0, 1, n_obs)
    d = (d_latent > np.median(d_latent)).astype(int)

    # Y: The continuous outcome. The true effect of D is 1.5.
    true_beta_d = 1.5
    y = 1.0 + true_beta_d * d + u

    # Create a DataFrame for analysis
    data = pd.DataFrame({'Y': y, 'D': d, 'Z1': z1, 'Z2': z2})
    data['Intercept'] = 1

    # Step 2: Perform the Hansen J-test
    # We estimate an IV model using BOTH Z1 and Z2 as instruments for D.
    # The model is overidentified because we have 2 instruments for 1 endogenous variable.
    # The Hansen J-test checks if the instruments are jointly valid.
    print("--- Testing Instrument Validity with Hansen J-Test ---")
    print("\nScenario: We are given instruments Z1 and Z2 for treatment D.")
    print("We know for sure that Z1 is a valid instrument.")
    print("We want to test if Z2 is also a valid instrument.")
    print("\nNull Hypothesis (H0) of the J-test: All instruments (Z1, Z2) are valid.")
    print("Alternative Hypothesis (H1): At least one instrument is invalid.")
    print("\nIf we reject H0, and given our assumption that Z1 is valid, we must conclude that Z2 is invalid.")

    model = IV2SLS(dependent=data.Y,
                   exog=data.Intercept,
                   endog=data.D,
                   instruments=data[['Z1', 'Z2']])

    results = model.fit(cov_type='robust')

    # Step 3: Output the test results and interpretation
    print("\n--- J-Test Results ---")
    # The 'sargan' attribute in the results object holds the test
    sargan_test = results.sargan
    print(sargan_test)

    # Automatically interpret the result based on a 5% significance level
    p_value = sargan_test.pval
    print("\n--- Interpretation ---")
    if p_value < 0.05:
        print(f"The p-value is {p_value:.4f}, which is less than 0.05.")
        print("We REJECT the null hypothesis. The test indicates that the instruments are not all valid.")
        print("Since we assumed Z1 was valid, this is strong evidence that Z2 is an INVALID instrument.")
    else:
        print(f"The p-value is {p_value:.4f}, which is not less than 0.05.")
        print("We FAIL to reject the null hypothesis. The test does not provide evidence against Z2's validity.")

    # Step 4: Output the estimated equation coefficients
    print("\n--- Estimated Final Equation (Y ~ 1 + D) ---")
    print("The numbers for the final equation are presented in the model summary below:")
    print(results)


# Execute the function
perform_instrument_validity_test()