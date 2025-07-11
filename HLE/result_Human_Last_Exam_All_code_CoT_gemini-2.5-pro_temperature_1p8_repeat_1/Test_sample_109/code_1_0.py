import numpy as np
import pandas as pd
from linearmodels.iv import IV2SLS
import statsmodels.api as sm

def run_instrument_validity_test():
    """
    This function simulates data for an IV setting where one instrument is valid
    and another is invalid. It then performs the Hansen J-test to test the
    joint validity of the instruments.
    """
    # 1. Set up simulation parameters
    n_obs = 1000
    np.random.seed(123)
    true_causal_effect = 2.0

    # 2. Generate data
    # 'c' is an unobserved confounder that affects both treatment D and outcome Y
    c = np.random.randn(n_obs)

    # Z1 is a valid instrument: it is not correlated with the confounder 'c'
    Z1 = np.random.randn(n_obs)

    # Z2 is an INvalid instrument: it is correlated with the confounder 'c'
    Z2 = 0.6 * c + np.random.randn(n_obs)

    # The error term 'u' in the outcome equation is driven by the confounder 'c'
    u = 1.2 * c + np.random.randn(n_obs)

    # The latent propensity for treatment depends on instruments and the confounder
    # This makes the treatment D endogenous because it's correlated with 'u' via 'c'
    d_latent = 0.7 * Z1 + 0.7 * Z2 - 0.8 * c + np.random.randn(n_obs)
    D = (d_latent > 0).astype(int)  # Binary treatment

    # The outcome Y is determined by the treatment D and the error 'u'
    Y = 1.5 + true_causal_effect * D + u

    # Create a pandas DataFrame
    data = pd.DataFrame({'Y': Y, 'D': D, 'Z1': Z1, 'Z2': Z2})
    data = sm.add_constant(data, prepend=False) # Add an intercept

    # 3. Define the IV model with both instruments (overidentified)
    dependent = data['Y']
    exog = data['const']      # Exogenous regressors (only intercept)
    endog = data['D']         # Endogenous regressor
    instruments = data[['Z1', 'Z2']] # Both valid and invalid instruments

    # 4. Fit the 2SLS model
    print("Fitting IV model using both Z1 (valid) and Z2 (invalid) as instruments...")
    iv_model = IV2SLS(dependent, exog, endog, instruments).fit(cov_type='robust')

    # 5. Perform and interpret the Hansen J-test (test of overidentifying restrictions)
    # The `sargan` attribute in linearmodels performs this test.
    # Null Hypothesis: All instruments are valid (exogenous).
    # We expect to REJECT this null hypothesis because Z2 is invalid.
    sargan_test = iv_model.sargan
    print("\n--- Hansen J-Test of Overidentifying Restrictions ---")
    print(f"Test Statistic: {sargan_test.stat:.4f}")
    print(f"P-value: {sargan_test.pval:.4f}")
    print(f"Degrees of Freedom: {sargan_test.df}")

    if sargan_test.pval < 0.05:
        print("\nConclusion: The p-value is less than 0.05.")
        print("We REJECT the null hypothesis that all instruments are valid.")
        print("This test correctly indicates that at least one of the instruments is invalid.")
    else:
        print("\nConclusion: The p-value is greater than 0.05.")
        print("We FAIL to reject the null hypothesis that all instruments are valid.")
        
    print("\nFor comparison, the estimated effect on D was:")
    print(iv_model.params)


if __name__ == '__main__':
    run_instrument_validity_test()
