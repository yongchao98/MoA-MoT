import numpy as np
import pandas as pd
from linearmodels.iv import IV2SLS

def run_validity_test_simulation():
    """
    Simulates data to demonstrate the Hansen J-test for instrument validity.

    In this simulation:
    - Y is the outcome, D is the endogenous treatment.
    - z1 is a valid instrument (correlated with D, uncorrelated with the error u).
    - z2 is an invalid instrument (correlated with D and the error u).
    - We test the null hypothesis that *both* z1 and z2 are valid.
    - A rejection of the null implies that at least one instrument is invalid.
    """
    # Set seed for reproducibility
    np.random.seed(42)
    
    # 1. Define simulation parameters
    n_obs = 1000  # Number of observations
    beta = 2.0    # True causal effect of D on Y

    # 2. Generate the data
    # Generate the structural error term
    u = np.random.normal(0, 1, n_obs)
    
    # Generate the instruments
    # z1 is a valid instrument, drawn from a standard normal distribution
    z1 = np.random.normal(0, 1, n_obs)
    
    # z2 is an INvalid instrument. We make it correlated with the error term 'u'.
    # This violates the exogeneity/exclusion restriction.
    z2 = 0.7 * u + np.random.normal(0, 1, n_obs)
    
    # Generate the endogenous treatment D
    # D depends on both instruments and is also correlated with u, making it endogenous.
    d_err = np.random.normal(0, 1, n_obs)
    D = 1.5 * z1 + 1.0 * z2 + 0.5 * u + d_err
    
    # Generate the final outcome Y
    # Y = beta*D + u
    Y = beta * D + u

    # 3. Prepare data for regression
    data = pd.DataFrame({
        'Y': Y,
        'D': D,
        'z1': z1,
        'z2': z2
    })
    # Add a constant for the intercept
    data['const'] = 1
    
    # 4. Fit the overidentified IV model and run the Hansen J-test
    # We use both z1 and z2 as instruments for the single endogenous variable D.
    # The model is overidentified (2 instruments, 1 endogenous regressor).
    model = IV2SLS(dependent=data['Y'],
                   exog=data['const'],
                   endog=data['D'],
                   instruments=data[['z1', 'z2']])
    
    results = model.fit(cov_type='robust')
    
    # 5. Print and interpret the results
    print("--- Over-identified IV Regression Results ---")
    print(results)
    
    j_stat = results.sargan.stat
    j_pvalue = results.sargan.pval
    
    print("\n--- Hansen J-Test for Overidentifying Restrictions ---")
    print(f"Hansen J-Statistic: {j_stat:.4f}")
    print(f"P-value: {j_pvalue:.4f}")
    
    print("\n--- Interpretation ---")
    print("The null hypothesis of the Hansen J-test is that all instruments are valid (exogenous).")
    if j_pvalue < 0.05:
        print("Result: The p-value is less than 0.05, so we REJECT the null hypothesis.")
        print("Conclusion: Assuming z1 is valid, this test provides strong evidence that z2 is an invalid instrument.")
    else:
        print("Result: The p-value is greater than 0.05, so we FAIL to reject the null hypothesis.")
        print("Conclusion: Assuming z1 is valid, there is not enough evidence to conclude that z2 is invalid.")

# Run the simulation and print the output
run_validity_test_simulation()