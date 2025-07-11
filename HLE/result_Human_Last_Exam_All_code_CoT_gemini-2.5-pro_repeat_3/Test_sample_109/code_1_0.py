import numpy as np
import pandas as pd
from linearmodels.iv import IV2SLS
import statsmodels.api as sm

def run_hansen_j_test_simulation():
    """
    Simulates data to demonstrate the Hansen J-test for instrument validity.

    In this simulation:
    - Y is the outcome.
    - D is the endogenous binary treatment.
    - Z1 is a valid instrument (uncorrelated with the error term U).
    - Z2 is an invalid instrument (correlated with the error term U).
    - The goal is to use the Hansen J-test to detect the invalidity of Z2,
      assuming Z1 is known to be valid.
    """
    # Set seed for reproducibility
    np.random.seed(42)
    n_obs = 1000

    # 1. Generate the data
    # Z1 is a valid instrument
    Z1 = np.random.normal(0, 1, n_obs)
    
    # Z2 will be an invalid instrument. We construct the error U to be correlated with Z2.
    Z2 = np.random.normal(0, 1, n_obs)
    
    # U is the error term in the outcome equation.
    # It is correlated with Z2, making Z2 an invalid instrument (violates exclusion).
    # It is uncorrelated with Z1, making Z1 a valid instrument.
    U = 0.8 * Z2 + np.random.normal(0, 1, n_obs)

    # 2. Create the endogenous treatment D
    # D depends on both instruments (relevance condition) and the error term U (endogeneity).
    # We use a latent variable approach to generate a binary D.
    d_latent = 0.5 * Z1 + 0.5 * Z2 + 0.6 * U + np.random.normal(0, 1, n_obs)
    D = (d_latent > 0).astype(int)

    # 3. Create the outcome Y
    # The true causal effect of D on Y is 2.0
    true_beta = 2.0
    Y = 1.0 + true_beta * D + U

    # 4. Prepare data for regression
    data = pd.DataFrame({'Y': Y, 'D': D, 'Z1': Z1, 'Z2': Z2})
    data['const'] = 1

    # 5. Perform IV (2SLS) regression
    # We use both Z1 and Z2 as instruments for the endogenous variable D.
    # The model is overidentified (2 instruments > 1 endogenous variable).
    model = IV2SLS(dependent=data['Y'],
                   exog=data['const'],
                   endog=data['D'],
                   instruments=data[['Z1', 'Z2']])
    
    results = model.fit(cov_type='robust')

    # 6. Print the results and explain the J-test
    print("--- IV (2SLS) Regression Results ---")
    print(results)

    j_stat = results.j_stat.stat
    j_pval = results.j_stat.pval

    print("\n--- Hansen J-Test of Overidentifying Restrictions ---")
    print(f"The Hansen J-statistic is: {j_stat:.4f}")
    print(f"The p-value of the J-statistic is: {j_pval:.4f}")

    print("\n--- Interpretation ---")
    print("The null hypothesis (H0) for the Hansen J-test is that all instruments are valid (exogenous).")
    print("In our setup, we maintain that Z1 is valid and test the validity of the additional instrument, Z2.")
    if j_pval < 0.05:
        print(f"Since the p-value ({j_pval:.4f}) is less than 0.05, we reject the null hypothesis.")
        print("This suggests that at least one instrument is invalid. Given our assumption that Z1 is valid,")
        print("we conclude that the instrument Z2 is invalid, which is correct by our simulation design.")
    else:
        print(f"Since the p-value ({j_pval:.4f}) is not less than 0.05, we fail to reject the null hypothesis.")
        print("This would suggest that we do not have enough evidence to claim that Z2 is invalid.")

# Run the simulation and print the output
run_hansen_j_test_simulation()