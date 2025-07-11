import numpy as np
import pandas as pd
from linearmodels.iv import IV2SLS

def perform_instrument_validity_test():
    """
    This function demonstrates the use of the Hansen J-test to check for
    instrument validity.

    We simulate a scenario with:
    - Y: outcome
    - D: endogenous treatment
    - Z1, Z3: valid instruments
    - Z2: an invalid instrument (correlated with the structural error)

    The Hansen J-test (or Sargan's test in this implementation) tests the null
    hypothesis that all instruments are exogenous (valid). A rejection of the null
    (low p-value) suggests at least one instrument is invalid.
    """
    # 1. Simulate Data
    np.random.seed(42)
    n_obs = 2000

    # Z1 and Z3 are valid instruments.
    # Z2 will be our INVALID instrument.
    z1 = np.random.normal(size=n_obs)
    z2 = np.random.normal(size=n_obs)
    z3 = np.random.normal(size=n_obs)

    # The structural error `u` is correlated with the invalid instrument `Z2`.
    # This violates the exogeneity condition for Z2.
    u = 0.6 * z2 + np.random.normal(size=n_obs)

    # The first-stage equation determines the endogenous variable D.
    # D must be correlated with all instruments to satisfy the relevance condition.
    d_latent = 0.5 + 1.0 * z1 + 1.0 * z2 + 1.0 * z3 + np.random.normal(size=n_obs)
    # Make D binary for this example
    d = (d_latent > np.mean(d_latent)).astype(int)

    # The structural (outcome) equation.
    y = 2.0 + 3.0 * d + u

    # Create a pandas DataFrame
    data = pd.DataFrame({
        'Y': y,
        'D': d,
        'Z1': z1,
        'Z2': z2,
        'Z3': z3
    })
    data['const'] = 1

    # 2. Perform Test 1: Using only VALID instruments (Z1, Z3)
    # The model is overidentified (2 instruments for 1 endogenous variable).
    # We expect the J-test to NOT reject the null of instrument validity.
    print("--- Test Case 1: Using only VALID instruments (Z1, Z3) ---")
    iv_valid = IV2SLS.from_formula('Y ~ 1 + [D ~ Z1 + Z3]', data=data).fit(cov_type='unadjusted')

    # The test is often called Sargan's J-test
    j_stat_valid = iv_valid.sargan.stat
    p_val_valid = iv_valid.sargan.pval

    print(f"Model: Y is outcome, D is endogenous, Z1 and Z3 are instruments.")
    print(f"Hansen J-statistic: {j_stat_valid:.4f}")
    print(f"P-value: {p_val_valid:.4f}")
    if p_val_valid > 0.05:
        print("Result: We FAIL to reject the null hypothesis. The test provides no evidence that the instruments are invalid.\n")
    else:
        print("Result: We REJECT the null hypothesis. The test suggests the instruments are invalid.\n")


    # 3. Perform Test 2: Using both VALID and INVALID instruments (Z1, Z2, Z3)
    # The model is overidentified (3 instruments for 1 endogenous variable).
    # We expect the J-test to REJECT the null because Z2 is invalid.
    print("--- Test Case 2: Using VALID (Z1, Z3) and INVALID (Z2) instruments ---")
    iv_invalid = IV2SLS.from_formula('Y ~ 1 + [D ~ Z1 + Z2 + Z3]', data=data).fit(cov_type='unadjusted')

    j_stat_invalid = iv_invalid.sargan.stat
    p_val_invalid = iv_invalid.sargan.pval

    print(f"Model: Y is outcome, D is endogenous, Z1, Z2, and Z3 are instruments.")
    print(f"Hansen J-statistic: {j_stat_invalid:.4f}")
    print(f"P-value: {p_val_invalid:.4f}")
    if p_val_invalid > 0.05:
        print("Result: We FAIL to reject the null hypothesis. The test provides no evidence that the instruments are invalid.\n")
    else:
        print("Result: We REJECT the null hypothesis. The test correctly indicates that at least one instrument in the set is invalid.\n")

if __name__ == '__main__':
    perform_instrument_validity_test()